# Local model runtime

L3's reasoning step runs an open-weights model **on your own machine** (decision
[D17](../mngmt/decisions.md)). Nothing is sent to a hosted service, which is why the step
refuses any endpoint that is not loopback — a local model behind a routable socket is a
hosted API, and nothing else about the call would look different.

This page is how to stand that up. Everything below is public software and public weights;
none of it touches the gated dataset.

The rest of the pipeline does **not** need any of this. L0–L2, the L3 ranking and L4 run
without a model. If the server is absent, L3 writes its ranking and records the gap in
`results/l3/integration.json` under `field_contract.reasoning`.

---

## What it needs

| | |
|---|---|
| Model | A 7–8B instruct model in GGUF, Q4_K_M |
| VRAM | ~5.4 GB free for full GPU offload at 8k context |
| Server | llama.cpp's `llama-server`, OpenAI-compatible endpoint |
| Bind | `127.0.0.1` only |

Measured on a GTX 1060 6 GB (driver 528.49) with Qwen2.5-7B-Instruct Q4_K_M, llama.cpp
b11065 Vulkan, full offload at 8192 context:

- prompt processing **181 tok/s** at 512, **175 tok/s** at 1536
- generation **24.0 tok/s** at 128, **23.7 tok/s** at 256
- peak **5,370 MiB** used, 672 MiB free
- three identical requests at `temperature 0`, `seed 42` returned byte-identical output

A pass over 83 candidates, two passes each, is roughly **1.6 hours**. This is a batch
pipeline, not an application.

---

## Pick a backend

**Use Vulkan unless you know you need CUDA.** It is the portable choice and, on older
cards, the only one that works.

Current llama.cpp CUDA builds no longer ship Pascal (compute 6.1) kernels: `llama-cli`
aborts with `invalid device function` inside `ggml_cuda_kernel_can_use_pdl`, a Hopper-era
code path. The CUDA 13.x builds fail for a blunter reason — CUDA 13 dropped Pascal
outright. Vulkan runs the same model at the same speed with no CUDA dependency.

If you use a CUDA build, note the runtime DLLs ship as a **separate** `cudart-*` download
in the same release. Without it the backend loads but reports no devices at all.

---

## Install

Download a release from [llama.cpp](https://github.com/ggml-org/llama.cpp/releases) and a
7–8B Q4_K_M GGUF. The spike used `Qwen/Qwen2.5-7B-Instruct-GGUF` (4.36 GB, split in two
files; point the server at part 1 and it finds the rest).

Do **not** fetch the model through `huggingface-cli` with this project's environment
loaded: `HF_HOME` points inside the custody root, and public model weights must not land
in a patient-data location. Download directly instead.

Verify the GPU is visible before downloading gigabytes of weights:

```
llama-cli --list-devices
```

You want a line naming your GPU. `(none)` means the backend did not load.

---

## Run

```
llama-server -m <model>.gguf -ngl 99 -c 8192 --host 127.0.0.1 --port 8080 --no-webui
```

`--host 127.0.0.1` is not optional. The pipeline refuses anything else, and the refusal is
the point: it is what makes "nothing leaves this machine" a property rather than a promise.

Check it:

```
curl -s http://127.0.0.1:8080/health
```

---

## If the pipeline runs in WSL and the server on Windows

WSL2's default NAT networking gives WSL its own network namespace, so `127.0.0.1` inside
WSL is **not** the Windows loopback. A server bound to Windows loopback is unreachable
from WSL, and so is the host's routed address.

Binding the server to `0.0.0.0` would "fix" this and is the wrong answer: it exposes the
model to anything that can route to the machine, and turns D17's guarantee into a claim.

Use **mirrored networking** instead, which makes WSL share the Windows network namespace
so that `127.0.0.1` means the same thing on both sides. Create `%USERPROFILE%\.wslconfig`:

```ini
[wsl2]
networkingMode=mirrored
```

Then `wsl --shutdown` and reopen your shell. This needs WSL 2.0.0+ and Windows 11 22H2+.
To undo it, delete the file and shut WSL down again.

Verify from inside WSL:

```
curl -s -o /dev/null -w '%{http_code}\n' http://127.0.0.1:8080/health   # expect 200
```

Building llama.cpp inside WSL is the other option, and is what a Linux-only judge would
do. The prebuilt Linux binaries will not help on Ubuntu 22.04: they are built against
Ubuntu 24.04 and need GLIBC 2.38, while 22.04 ships 2.35.

---

## Configure the pipeline

In `config/pipeline.yaml`:

```yaml
model: qwen2.5-7b-instruct-q4_k_m      # recorded in results/_manifest.json
reasoning_endpoint: http://127.0.0.1:8080

l3:
  reasoning: true
  reasoning_top_n: 100                 # scope before L4 has run
  reasoning_timeout_seconds: 300
```

The model id is free text and is recorded with every rationale — a rationale whose model
is unknown is not reproducible evidence, so set it to whatever you actually served.

Then:

```
python -m src.pipeline --only l3_integrate
```

Scope: after L4 has run, the step reasons over `results/l4/survivors.tsv` only, which is
what the dossier needs. Before that it takes the top `reasoning_top_n` of L3's ranking.
Which one was used is recorded in `results/l3/reasoning/reasoning.json`.

---

## What the step will and will not do

It asks the model only what a model of this size can answer from supplied text: whether
the evidence supports the mechanism, whether there is a contradiction in it, and which of
the supplied citation keys are load-bearing.

It does **not** ask the model to grade evidence. During the 2026-09-20 spike a 7B graded
cell-line work as `clinical`. The grade comes from channel E, which already graded its
literature, and travels with the claim.

It does **not** ask the model for a confidence number. Confidence is computed in Python
from the supplied grade, the number of supporting channels, whether that convergence is
discriminating, the adverse-direction records channel E counted, and the model's rubric.
That makes it testable. It is a hypothesis-strength score and never a probability of
clinical benefit.

Every citation key the model returns is checked against the packet and dropped if it was
not supplied. Drops are recorded in `rationales.json`, not hidden: a model inventing
citations is a finding about the model.

# PICurv Storage Package Guide

Protect, offload, restore, and reclaim PICurv artifacts through an `rclone` remote.
The package is split so that each file answers one question, and the dependency
direction runs strictly down this list — nothing here imports from a module below it.

| Module | Question it answers |
|---|---|
| `models.py` | What are the constants, the profile, and the local state marker? |
| `transport.py` | How does a byte reach or leave the remote? |
| `safety.py` | What makes an artifact unsafe to touch right now? |
| `inventory.py` | What does this artifact contain, and which target was selected? |
| `packaging.py` | How is that content turned into chunks, and back? |
| `compatibility.py` | What provenance is captured, and how do paths survive a restore? |
| `catalog.py` | What is on the remote, and what still references it? |
| `operations.py` | Plan, protect, offload, restore, prune, and the command surface. |

`__init__.py` re-exports the whole surface, including the private helpers the
conductor and the test suite address, so importers write `picurv_cli.storage.X`
without knowing which file `X` lives in.

Two functions are reached through their module rather than imported by name —
`transport._run_rclone` and `transport._read_remote_bytes`. They are the process
boundary, and routing every caller through the module means a test replaces the
boundary in one place instead of once per importer.

## Where to add something

- A new refusal belongs in `safety.py`, not inline at the call site, so that
  `plan` and the real operation cannot disagree about what is safe.
- A new thing to inventory belongs in `inventory.py` and needs a component in
  `_classify_component`. An unclassified path is archived and never pruned, which
  is deliberate; leaving something unclassified is a choice, not an oversight.
- A new remote layout belongs behind `transport._chunk_remote_path`, which is the
  one place that decides where a chunk's payload lives.

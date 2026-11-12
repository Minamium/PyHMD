# Agent skills

[English](skills.md) | [日本語](ja/skills.md) | [Index](index.md)

Two repository-scoped skills live in `.agents/skills/`. Keeping them in the
checkout makes their guidance follow the same version as the code. Skill bodies
and UI metadata are in English; replies should follow the user's language.

| Skill | Use it for | Entry point |
| --- | --- | --- |
| `$pyhmd-develop` | Source changes, numerical/API contracts, tests, packaging, documentation | [SKILL.md](../.agents/skills/pyhmd-develop/SKILL.md) |
| `$pyhmd-use` | Job status, time estimates, settings, runs, restart, diagnostics, plots | [SKILL.md](../.agents/skills/pyhmd-use/SKILL.md) |

In a chat working in this checkout, invoke the skill by name. If the client has
not loaded repository skills, supply the linked `SKILL.md` path explicitly.
No global skill installation or unrelated user configuration is modified.

## Example requests

```text
Use $pyhmd-develop to change the result manifest API and test restart compatibility.
Use $pyhmd-develop to refactor the CUDA reconstruction kernel without changing the numerics.
Use $pyhmd-use to prepare a 600×128 thin-disk run on an existing single-GPU allocation.
Use $pyhmd-use to check this PBS job and estimate its time to t=500 from comparable diagnostics.
Use $pyhmd-use to inspect this completed run and plot the inner disk and transport diagnostics.
```

The development skill routes changes to implementation owners and appropriate
checks. It separates refactors from deliberate numerical changes and preserves
bundled-source, output, and checkpoint contracts. The usage skill operates the
existing interfaces and distinguishes completed, partial, failed, numerically
validated, and scientifically reproduced results.

Neither skill contains private cluster connections or grants extra authority to
allocate compute, cancel unrelated jobs, or publish code. Existing task-specific
authorization still applies; an invocation is not an extra approval requirement.

## Maintenance

Read only the task-relevant references linked by the skill. Add supporting files
only when they provide reusable information not already maintained in the main
documentation. Keep the two discovery descriptions distinct. Preserve default
automatic selection; explicit invocation remains available.

Run `python3 scripts/check_docs.py` for local links and language pairs. When
editing a skill, also run the skill-creator `quick_validate.py` helper available
in the authoring environment and review `agents/openai.yaml`. That helper is not
a runtime dependency of PyHMD. Do not claim its syntax checks prove behavior.

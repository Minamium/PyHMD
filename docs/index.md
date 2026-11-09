# Documentation

[English](index.md) | [日本語](ja/index.md) | [README](../README.md)

English pages are primary. Japanese counterparts retain the same supported
features, commands, and scientific limitations. Identifiers and machine-readable
artifacts are not translated.

| Guide | Read it for | 日本語 |
| --- | --- | --- |
| [Architecture](architecture.md) | Ownership and module boundaries | [設計](ja/architecture.md) |
| [kineticEQ design](kineticeq-design.md) | Concrete design decisions inherited from kineticEQ | [設計の対応](ja/kineticeq-design.md) |
| [Configuration](configuration.md) | Typed settings, JSON, INI, preset merging | [設定](ja/configuration.md) |
| [Results](results.md) | Artifacts, terminal states, strict restart | [結果と再開](ja/results.md) |
| [HPC](hpc.md) | Build targets, MPI launch, PBS/Slurm native execution | [HPC実行](ja/hpc.md) |
| [Scientific status](scientific-status.md) | Mathematical scope and known reproduction gaps | [科学的検証](ja/scientific-status.md) |
| [Validation](validation.md) | Test commands, evidence, hardware limits | [検証](ja/validation.md) |
| [Agent skills](skills.md) | Development and simulation workflows | [skills](ja/skills.md) |

The detailed [low-β thermal regression notes](../tests/thermal_robustness.md)
remain an English technical test record. Machine-readable evidence lives in
[reference_import.json](reference_import.json) and
[migration_cpu.json](validation/migration_cpu.json); these are shared by both languages.

## Maintaining the two languages

Use `README.md` and `docs/*.md` for English; use `README.ja.md` and
`docs/ja/*.md` for Japanese. Update behavior, examples, and limitations in both
versions together. Link counterpart pages and keep command/configuration names
identical. Run `python3 scripts/check_docs.py` before committing documentation.

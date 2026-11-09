# ドキュメント

[English](../index.md) | [日本語](index.md) | [README](../../README.ja.md)

英語を主版とし、日本語版でも対応機能・コマンド・科学的な制限を揃えます。
API名、設定キー、機械可読ファイルの識別子は英語のままです。

| ガイド | 内容 | English |
| --- | --- | --- |
| [設計](architecture.md) | モジュールと所有権 | [Architecture](../architecture.md) |
| [kineticEQとの対応](kineticeq-design.md) | 参考にした具体的な設計 | [Design mapping](../kineticeq-design.md) |
| [設定](configuration.md) | 型、JSON、INI、preset | [Configuration](../configuration.md) |
| [結果と再開](results.md) | 保存物、終了状態、再開条件 | [Results](../results.md) |
| [HPC](hpc.md) | ビルドとMPI、PBS/Slurmからの実行 | [HPC](../hpc.md) |
| [科学的な検証範囲](scientific-status.md) | 計算モデルと再現の未確認事項 | [Scientific status](../scientific-status.md) |
| [検証](validation.md) | テスト、実測済み範囲、制限 | [Validation](../validation.md) |
| [開発用・利用用skill](skills.md) | エージェント向けの作業手順 | [Agent skills](../skills.md) |

[低β熱回復テストの詳細記録](../../tests/thermal_robustness.md)は英語の技術資料です。
[reference_import.json](../reference_import.json)と
[migration_cpu.json](../validation/migration_cpu.json)は両言語で共有する機械可読の検証記録です。

英語は`README.md`と`docs/*.md`、日本語は`README.ja.md`と`docs/ja/*.md`で管理します。
挙動・コマンド例・制限事項を変えるときは両版を更新し、相互リンクを付けます。
確認コマンドは`python3 scripts/check_docs.py`です。

- [Phase 0：計測・検証基盤](phase0.md)（[English](../phase0.md)）：計画の検証、probe、凍結対照と未完了事項。

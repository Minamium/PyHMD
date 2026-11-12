# 開発用・利用用skill

[English](../skills.md) | [日本語](skills.md) | [ドキュメント一覧](index.md)

2つのskillを`.agents/skills/`に同梱しています。コードと同じ版で管理するため、
このcheckoutを作業対象として使います。skill本文とUI情報は英語ですが、返信言語はユーザに合わせます。

| 名前 | 担当 | 本文 |
| --- | --- | --- |
| `$pyhmd-develop` | 実装、数値/APIの契約、テスト、配布、文書 | [SKILL.md](../../.agents/skills/pyhmd-develop/SKILL.md) |
| `$pyhmd-use` | ジョブ状態、所要時間、設定、実行、再開、診断、可視化 | [SKILL.md](../../.agents/skills/pyhmd-use/SKILL.md) |

このcheckoutのチャットで名前を指定します。クライアントがリポジトリのskillをまだ読み込んでいなければ、
上記`SKILL.md`のパスを明示してください。グローバルなskill配置やユーザ設定は変更しません。

```text
$pyhmd-develop を使って結果manifestのAPIを修正し、再開の互換性を検証して。
$pyhmd-develop を使ってCUDAの再構成kernelを数値式を変えずに整理して。
$pyhmd-use を使って、確保済み単GPUで600×128の円盤計算を実行する設定を用意して。
$pyhmd-use を使って、このPBSジョブの状態を確認し、同等条件の診断からt=500までの所要時間を見積もって。
$pyhmd-use を使って、この結果の内部円盤と輸送診断を可視化して。
```

開発用は変更の責務と検証方法を案内し、構造変更と数値変更を分けます。
利用用は既存APIを使い、完走・途中終了・失敗・数値整合性・科学的再現を区別します。
どちらも接続情報を含まず、追加資源の確保、無関係なジョブ停止、公開の権限を新しく与えません。
現在の依頼で得ている権限はそのまま適用され、skillを呼ぶだけで追加確認が必要になることもありません。

変更時は関連する参照だけを読み、既存文書を複製せず、2つの用途と選択条件を区別します。
自動選択の既定値を維持しつつ、明示指定も使えます。
`python3 scripts/check_docs.py`に加え、作成環境のskill-creatorにある`quick_validate.py`で
構文を検証し、`agents/openai.yaml`も確認してください。このhelperはPyHMDの実行依存ではなく、
構文検証だけでskillの実際の判断を保証するものでもありません。

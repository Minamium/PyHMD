# Numerical and scientific scope

PyHMDは現在、SF2014を対象に開発している研究用の再実装です。
構造整理で既存の物理設定・対応座標・数値演算を増減させることは意図していません。

対応する式は古典圧縮性抵抗性MHD、理想気体、球面の点質量重力です。
PPM/HLL/SSP-RK3と、軸対称では磁束関数Ψを使うCT、異方性Ohmic項を陽的に更新します。
全体の精度次数は、境界・抵抗項・熱回復を含めて独立に評価する必要があります。
PPM/RK3という名前だけで全問題の三次収束は主張しません。

## Known fidelity gaps

- 元のPLUTO4.0のinit.c/definitions.h/pluto.ini/res_eta.cは未入手です。
- 内側のセル中心磁場角と現在のΨ境界処方は同一の離散条件ではありません。
- ghost速度の符号と実HLL質量流束の向きは別であり、sinkの物理的意図の照合が必要です。
- 初期ガス平衡のhydrostatic/magnetohydrostaticの解釈差があります。
- compatible Ohmic cooling、CT補正、内部エネルギー回復の原著との同一性は未確認です。
- 標準拡散モデルとstrongモデルは別実験です。長時間安定モデルへ自動切替しません。
- CPU検証をCUDA検証へ、単一ノードMPIの結果を複数ノードへ一般化しません。

完走、divBの小ささ、補正込みの収支閉合だけでは物理解を保証できません。
この版の`completed`と数値validator合格は、SF2014の再現成功を意味しません。

参考：[Stepanovs & Fendt 2014](https://arxiv.org/html/1407.3136)、
[metric PPMの幾何](https://arxiv.org/abs/1404.0537)、
[UCTの構成](https://doi.org/10.1016/j.jcp.2020.109748)。

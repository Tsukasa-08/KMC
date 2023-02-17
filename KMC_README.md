# KMC
拡散係数を計算するためのKMCプログラム

## 入力ファイル

<details>
<summary>INPUT</summary>

- KMCの設定ファイル
- VASPのINCARのように、”[設定したい変数]=[設定値]”と記述する

```
#設定ファイル

#Monte Carlo Step per Particle (MCSP, 修論:Nsteps) : 1度のKMCシミュレーションで1粒子あたり平均何回ジャンプをさせるか
MCSP=10000

#統計平均を取る粒子数 (修論:Nparticles)
AVERAGE=1

#系内に配置する拡散粒子数
NDIFFS=108

#電場を設定するか否か (電場なし=0, 電場あり=1)
EFIELDON=0

#電場によるジャンプ頻度補正の目安、10の対数値で入力する (修論 : アルファの指数に対応)
CORRECT=-1.6

#電場をかける軸方向(+x=1, +y=2, +z=3, -x=-1, -y=-2, -z=-3)
AXIS=1

#PES勾配を打ち消す場合の電場ベクトル
ANTIDRIFT = [0.0, 0.0, 0.0]

#ジャンプリストの平均ジャンプ距離(Ang.) いずれはプログラムで抽出できた方がいいかも
DISTANCEJUMP=1.5

#シミュレーション温度(K)
TEMP=600

#sitePE.datを読み込むかどうか(yes=1, no=0)
SITEPEREAD=1

#blocking_list.csvを読み込むかどうか(yes=1, no=0)
BLOCKINGLISTREAD=1

#blocking_listによるarea blockingを行うかどうか(yes=1, no=0)
BLOCKING=1

#hoppingやrotationのcounterを出力するか(BaZrO3専用)
ROTHOPCOUNT=0

#拡散を考える次元(default=3),2次元PESでも3で問題ない?
DIMENSIONALITY=3
```
</details>

<details>
<summary>JMPDATA</summary>

- 隣接するサイト間を結ぶジャンプ頻度を記述したファイル
- 本研究では、空孔を1、拡散粒子(=プロトン)を2として区別していた
- makejmpdata.shで豊浦先生および藤井さんのjmpdata.csvをJMPDATAの形式に変換可能

```
#InitialSiteID, InitialSiteAtom, FinalSiteID, FinalSiteID, Frequency[Hz}
1,2,973,1,21282711265566.7021
1,2,8281,1,21282711265566.7021
```
</details>

<details>
<summary>POSCAR</summary>

- 対象とする格子を記述したファイル、VASPに用いるファイルと同一のフォーマット
</details>

<details>
<summary>sitePE.dat</summary>

- 各サイトのエネルギーを記載したファイル
- POSCARの行と対応している必要あり
</details>

<details>
<summary>blocking_list.csv</summary>

- area blockingを考慮するサイトを1行ずつカンマ区切りで並べたファイル
</details>

## 出力ファイル

<details>
<summary>OUTPUT</summary>

- toml形式で作成されたファイル
- 計算条件などが利用しやすい形でまとまっている

```
#This in OUTPUT file written by toml format.

#total_time [s] (average of all KMCs : 全KMCシミュレーションの平均値を計算している, ステップ数を決めればほとんど一致するので)
total_time = 1.204e-08

#concentration [/Ang.^3] 体積あたりの拡散粒子数、濃度
concentration = 1.9487e-03

#temperture [K]
temperture = 600

#ion_charge 拡散種の電荷(integer)
ion_charge = 1

#mean_displacement [Ang.]
mean_displacement = [6.1, -1.1, 2.5]
```
</details>

<details>
<summary>log_cout</summary>

- 出力ログ
- 入力ファイルから読み取った情報、およびプログラムの経過や実行時間などが記録される
</details>

<details>
<summary>mean_displacement.csv</summary>

- 平均変位を出力するファイル

```
#the number of KMC, diffusion_id, displacement in x direction [Ang.], displacement in y direction [Ang.], displacement in z direction [Ang.], sum of squared displacement of each jumps in x, y, z [Ang.^2], start_SiteID, end_SiteID, jump_counter of total, rotation, hopping [times]
KMC_times,diffusion_id,dx,dy,dz,sum_x2,sum_y2,sum_z2,start_site,end_site,jump_counter,rot_counter,hop_counter
```
</details>


## 使用方法

1. KMC/base_dirを任意のディレクトリにダウンロードする
2. Makefileのあるディレクトリ内で、「make all」コマンドを実行しコンパイルを完了させる
3. 入力ファイルが用意されたディレクトリを用意する
4. 入力ファイル内で、C++の実行ファイルである「test_kmc_test」を実行 (プログラム自体はどこにあるものを呼び出してもよい)
5. ログや経過は標準出力に出力されるので、適宜 >& log_coutなどでファイルにリダイレクトする

## プログラムの修正を行う場合
- base_dir/src/kmc_test.cpp がソースコード
- base_dir/inc に必要なインクルードファイルはすべて入っている
- base_dir/src/Makefile に従ってmakeすればコンパイルされる


## 計算前の確認事項

- MCSPは所定の設定か
- POSCARとJMPDATAは対応しているか
- POSCARとblocking_list.csvは対応しているか
- POSCARとsitePE.datは対応しているか
- make_calc_directory.shによるsedの際にdefaultになるべき値はその通りか
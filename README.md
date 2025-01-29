# self-dual-PotBKZ

本ソースたちは北九州で開催のSCIS2025にて発表した研究成果のソースです．

論文題目は『自己双対型PotBKZ基底簡約の提案とBKZとの比較』，発表ブロックは[2D2](https://www.iwsec.org/scis/2025/program.html#2D2)です．

本リポジトリには，停止性の保証されたBKZアルゴリズム[SE94]の変種，PotBKZとその双対版である双対型PotBKZ，更に自己双対型となる自己双対型PotBKZが入っています．<br>具体的には，BKZ[SE94]はLLLアルゴリズム[LLL82]などと比較して非常に良い基底を出力することが実験的にも知られているが，トレードオフとして，格子次元に関して指数的で，停止性は保証されていません．<br>本リポジトリに含まれているPotBKZや自己双対型PotBKZは基底に関して定まる**ポテンシャル**という量を単調減少させることで，格子次元に関して多項式回のツアーで停止するBKZの変種です．

The sources in this repository are research results that was presented in SCIS2025 at Kita-kyushu city in Japan. The research results are new variants of BKZ with provably termination **PotBKZ** and its dual version and self-dual version.<br>Specifically, BKZ algorithm[SE94] outputs much better lattice basis than, for example, LLL algorithm[LLL82]. But as trade-off, BKZ algorithm is exponentially on lattice dimension and is not guaranteed terminating.<br>PotBKZ, its dual version, and its self-dual version are new variants of BKZ algorithm that terminates polynomaial tour on lattice dimension by monotonically decreasing potential of lattice basis.

- [LLL82] Arjen Klaas Lenstra and Hendrik Willem Lenstra and László Lovász, "Factoring polynomials with rational coefficients", 1982
- [SE94] Claus-Peter Schnorr and Martin Euchner, "Lattice basis reduction: Improved practical algorithms and solving subset sum problems", 1994

# self-dual-PotBKZ

## How to Install This
First, you need to clone this repository to your computer:
```shell
$ git clone https://github.com/satoshin-des/self-dual-PotBKZ.git
```

Using ``cd`` command to change the directory ``self-dual-PotBKZ``:
```shell
$ cd self-dual-PotBKZ
```


The sources in this repository are research results that was presented in SCIS2025 at Kita-kyushu city in Japan. The research results are new variants of BKZ with provably termination **PotBKZ** and its dual version and self-dual version.<br>Specifically, BKZ algorithm[SE94] outputs much better lattice basis than, for example, LLL algorithm[LLL82]. But as trade-off, BKZ algorithm is exponentially on lattice dimension and is not guaranteed terminating.<br>PotBKZ, its dual version, and its self-dual version are new variants of BKZ algorithm that terminates polynomaial tour on lattice dimension by monotonically decreasing potential of lattice basis.

- [LLL82] Arjen Klaas Lenstra and Hendrik Willem Lenstra and László Lovász, "Factoring polynomials with rational coefficients", 1982
- [SE94] Claus-Peter Schnorr and Martin Euchner, "Lattice basis reduction: Improved practical algorithms and solving subset sum problems", 1994

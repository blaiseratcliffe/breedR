# breedR on Linux

Two things get in the way of the BLUPF90 programs on Linux: the download from UGA fails at the moment, and most binaries need the Intel MKL runtime to start.

## Downloads from UGA fail (TLS)

Downloads from nce.ads.uga.edu fail on Linux with an SSL certificate error. This affects the download during package installation, and `install_progsf90()`, `install_genomic_programs()` and `install_renumf90()`. Windows is not affected. macOS is untested. Progress is tracked in [#116](https://github.com/blaiseratcliffe/breedR/issues/116).

The cause is the server's certificate chain, as checked on 2026-10-02:

```
nce.ads.uga.edu                                  (expires 2027-04-15)
  issued by: InCommon Intermediate CA - OVG2C
  issued by: emSign Root TLS CA - G1
```

The server sends only its own certificate, without the intermediate. Fetching the intermediate does not help: its issuer, `emSign Root TLS CA - G1`, is not in Mozilla's root store, so it is not in the `ca-certificates` bundle of Debian, Ubuntu and most other distributions, and the chain cannot be verified. You can check what the server sends with:

```sh
openssl s_client -connect nce.ads.uga.edu:443 -servername nce.ads.uga.edu -showcerts </dev/null
```

While this holds there is no safe general fix on the user's side. Adding `emSign Root TLS CA - G1` to the system trust store would make the download work, but it makes every program on the machine trust a root that Mozilla does not include. That is a security decision for you or your system administrator, not a step this page recommends.

The previous advice to trust `InCommon RSA Server CA 2` was for the certificate UGA used before its renewal, and no longer applies.

To install breedR without trying the download, set `BREEDR_SKIP_INSTALL_BINARIES=true` in the environment of `R CMD INSTALL`. The binaries can then be installed later with `install_progsf90()`, once the download works again.

## Intel MKL and OpenMP libraries

Most of the BLUPF90+ binaries (all but renumf90 and postgibbsf90) are linked against Intel MKL and the Intel OpenMP runtime, and will not start without them:

```
error while loading shared libraries: libmkl_intel_lp64.so.2
```

Intel's PyPI wheels provide the libraries (`libmkl_intel_lp64.so.2`, `libmkl_intel_thread.so.2`, `libmkl_core.so.2` and `libiomp5.so`):

```sh
python3 -m pip download --no-deps --only-binary=:all: --dest wheels \
  mkl==2024.2.2 intel-openmp==2024.2.1
mkdir -p ~/mkl && for w in wheels/*.whl; do
  unzip -q -j -o "$w" '*.data/data/lib/*' -d ~/mkl
done
export LD_LIBRARY_PATH=~/mkl${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}   # before starting R
```

Back to the [README](../README.md).

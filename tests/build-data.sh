#!/usr/bin/bash

set -euo pipefail
set -x

# Call from project root to re-build test data.
#
# $ bash tests/build-data.sh

mkdir -p tests/genome

cargo run genome --out-file tests/genome/expected.fa -l 100 -l 100
cargo run genome --out-file tests/genome/expected-ll-80.fa -l 100 -l 100 --line-length 80

mkdir -p tests/methylation

cp tests/genome/expected.fa tests/methylation/input.fa
cargo run methylation --in tests/methylation/input.fa --out tests/methylation/expected.fa

mkdir -p tests/frag-sequencing

cargo run genome --out-file tests/frag-sequencing/input-1.fa -l 10000
cargo run frag-sequencing --seq-technology sanger --in tests/frag-sequencing/input-1.fa --out tests/frag-sequencing/output-1-sanger.fa
cargo run frag-sequencing --seq-technology sanger --in tests/frag-sequencing/input-1.fa --out tests/frag-sequencing/output-1-sanger-read-info.fa --embed-read-info

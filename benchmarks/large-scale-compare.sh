#!/bin/sh
set -e

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work_root=$STM_LARGE_WORK_ROOT
if [ -z "$work_root" ]; then
  work_root=$(mktemp -d /tmp/stm-large-scale.XXXXXX)
fi
dependency_lib=$STM_DEPENDENCY_LIBRARY
makevars_user=$R_MAKEVARS_USER
seeds=$STM_BENCHMARK_SEEDS
if [ -z "$seeds" ]; then seeds=1; fi
documents=$STM_BENCHMARK_DOCUMENTS
if [ -z "$documents" ]; then documents=50000; fi
epochs=$STM_BENCHMARK_EPOCHS
if [ -z "$epochs" ]; then epochs=20; fi
corpus=$STM_BENCHMARK_CORPUS
if [ -z "$corpus" ]; then corpus=$work_root/anchored-$documents.rds; fi
results=$STM_BENCHMARK_RESULTS
if [ -z "$results" ]; then results=$work_root/results; fi

mkdir -p "$work_root/svi-source" "$work_root/current-library" \
  "$work_root/svi-library" "$results"

if [ ! -f "$corpus" ]; then
  (
    cd "$repo_root"
    Rscript benchmarks/large-scale-generate.R \
      "--documents=$documents" "--output=$corpus"
  )
fi

(
  cd "$repo_root"
  git archive origin/svi | tar -x -C "$work_root/svi-source"
)

install_package() {
  source_path=$1
  library_path=$2
  if [ -n "$dependency_lib" ]; then
    R_LIBS="$dependency_lib" R_MAKEVARS_USER="$makevars_user" \
      R CMD INSTALL --preclean --library="$library_path" "$source_path"
  else
    R_MAKEVARS_USER="$makevars_user" \
      R CMD INSTALL --preclean --library="$library_path" "$source_path"
  fi
}

install_package "$repo_root" "$work_root/current-library"
install_package "$work_root/svi-source" "$work_root/svi-library"

run_method() {
  method=$1
  seed=$2
  library_path=$3
  output=$results/$method-seed$seed.rds
  libs=$library_path
  if [ -n "$dependency_lib" ]; then
    libs=$libs:$dependency_lib
  fi
  (
    cd "$repo_root"
    R_LIBS="$libs" /usr/bin/time -p Rscript benchmarks/large-scale-run.R \
      "--method=$method" "--corpus=$corpus" "--output=$output" \
      "--seed=$seed" "--max_epochs=$epochs" \
      2>"$results/$method-seed$seed.time.txt"
  )
}

for seed in $seeds; do
  run_method batch "$seed" "$work_root/current-library"
  run_method memoized "$seed" "$work_root/current-library"
  run_method svi "$seed" "$work_root/svi-library"
done

(
  cd "$repo_root"
  Rscript benchmarks/large-scale-summarize.R \
    "--results=$results" "--output=$results/summary.csv"
)
printf 'Large-scale comparison saved under %s\n' "$work_root"

#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Answers "which modules do I need" empirically, instead of guessing.
# Safe, read-only, ~10 seconds. Run it on the login node, NOT via bsub:
#
#   cd ~/docker-testing && bash wexac_probe.sh 2>&1 | tee probe.txt
#
# Then send me probe.txt (or just paste the summary at the end).
# ---------------------------------------------------------------------------
echo "=============================================================="
echo " WEXAC PROBE   $(date)   $(hostname)   $(id -un)"
echo "=============================================================="

hdr () { echo; echo "--- $* ---"; }

hdr "1. module system"
if command -v module >/dev/null 2>&1; then
  echo "module command: present"
  module --version 2>&1 | head -3
else
  echo "module command: NOT on PATH"
  for f in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
    [ -f "$f" ] && echo "  but found init script: $f  (source it first)"
  done
fi

hdr "2. modules that matter here"
for pat in -i singularity apptainer python conda miniconda mamba R/ r- java meme; do
  :
done
if command -v module >/dev/null 2>&1; then
  module avail 2>&1 | grep -iE 'singular|apptainer|conda|mamba|^ *R/|python|meme|bowtie|samtools' | sed 's/^/  /' | head -40
  echo
  echo "  (full list: module avail 2>&1 | less)"
else
  echo "  skipped, no module command"
fi

hdr "3. currently loaded"
( module list ) 2>&1 | head -20

hdr "4. container tooling"
for t in singularity apptainer docker podman; do
  if command -v "$t" >/dev/null 2>&1; then
    printf '  %-12s YES  %s\n' "$t" "$(command -v "$t")"
    "$t" --version 2>&1 | head -1 | sed 's/^/               /'
  else
    printf '  %-12s no (on this node, before any module load)\n' "$t"
  fi
done
echo "  can you reach the build server?"
timeout 8 ssh -o BatchMode=yes -o StrictHostKeyChecking=no \
    docker-gpu.wexac.weizmann.ac.il 'echo "    docker-gpu reachable; docker: $(command -v docker || echo MISSING)"' \
    2>&1 | sed 's/^/    /' || echo "    could not ssh non-interactively (may just need a password/key)"

hdr "5. conda / python already available"
echo "  CONDA_PREFIX = ${CONDA_PREFIX:-<unset>}"
for t in conda mamba micromamba python3 Rscript; do
  if command -v "$t" >/dev/null 2>&1; then
    printf '  %-12s %s\n' "$t" "$(command -v "$t")"
  else
    printf '  %-12s no\n' "$t"
  fi
done
command -v python3 >/dev/null && python3 -c "import sys;print('  python', sys.version.split()[0])"
command -v Rscript  >/dev/null && Rscript -e 'cat("  R", R.version.string, "\n")' 2>/dev/null

hdr "6. LSF queues"
if command -v bqueues >/dev/null 2>&1; then
  bqueues 2>/dev/null | head -25 | sed 's/^/  /'
  echo
  echo "  docker application profiles:"
  bapp 2>/dev/null | sed 's/^/    /' | head -15 || echo "    bapp not available"
else
  echo "  bqueues not found - are you on a login node?"
fi

hdr "7. storage"
df -h . 2>/dev/null | sed 's/^/  /'
echo "  quota:"
( lfs quota -u "$(id -un)" . 2>/dev/null || quota -s 2>/dev/null || echo "    n/a" ) | sed 's/^/  /' | head -8

hdr "8. outbound network (needed to download from ENA)"
if timeout 15 curl -sfI https://www.ebi.ac.uk/ena/portal/api/ >/dev/null 2>&1; then
  echo "  ENA API reachable from this node"
else
  echo "  ENA API NOT reachable from the login node"
  echo "  -> downloads may need a proxy, or must run from a node with egress"
fi
timeout 15 curl -sfI https://ftp.sra.ebi.ac.uk >/dev/null 2>&1 \
  && echo "  ENA FTP/HTTPS reachable" || echo "  ENA FTP/HTTPS NOT reachable"

echo
echo "=============================================================="
echo " SUMMARY"
echo "=============================================================="
printf ' module system : %s\n' "$(command -v module >/dev/null && echo yes || echo 'no - source an init script')"
printf ' singularity   : %s\n' "$(command -v singularity >/dev/null && echo 'on PATH already' || echo 'needs: ml Singularity')"
printf ' conda         : %s\n' "${CONDA_PREFIX:-not active}"
printf ' python3       : %s\n' "$(command -v python3 || echo missing)"
printf ' Rscript       : %s\n' "$(command -v Rscript  || echo missing)"
printf ' bsub          : %s\n' "$(command -v bsub     || echo missing)"
echo
echo "Next: bash setup_fast.sh   (builds a conda env and runs the pipeline"
echo "      without waiting for a Docker image)"

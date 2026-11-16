# Building and running the image on Wexac

Source: <https://hpcwiki.weizmann.ac.il/en/home/lsf/container>

You cannot build containers on the Wexac compute or login nodes, building
needs root. Weizmann solves this with a **dedicated build server**:

```
docker-gpu.wexac.weizmann.ac.il
```

You have Docker privileges there and nowhere else. Build there, then either
push to the internal registry or copy a `.sif` to your home.

---

## Route A, Docker via LSF (default, `CONTAINER_MODE=docker-lsf`)

LSF launches the whole job inside the image using its `docker-cpu` application
profile. Nothing to convert, nothing to copy.

**1. Build and push, on the build server:**

```bash
ssh docker-gpu.wexac.weizmann.ac.il
cd /path/to/phip_pipeline
docker build -f docker/Dockerfile -t phip-pipeline:0.1 docker/
docker tag phip-pipeline:0.1 ops:5000/phip-pipeline:0.1
docker push ops:5000/phip-pipeline:0.1
```

Check it landed: `ops_repositories`

**2. Run from a Wexac access node**, this is what `run_all.sh` does:

```bash
bsub -q short -app docker-cpu -env LSB_CONTAINER_IMAGE=ops:5000/phip-pipeline:0.1 \
     -R rusage[mem=16000] -n 8 bash scripts/02_count.lsf
```

Interactive, to poke around inside:

```bash
bsub -q interactive -app docker-cpu-interactive \
     -env LSB_CONTAINER_IMAGE=ops:5000/phip-pipeline:0.1 -Is /bin/bash
```

Confirm you really are inside the image, `/.dockerenv` exists only in a
container:

```bash
[ -f /.dockerenv ] && echo "in container" || echo "on the host"
```

Set in `config/params.sh`:

```sh
CONTAINER_MODE="docker-lsf"
IMAGE="ops:5000/phip-pipeline:0.1"
```

---

## Route B, Singularity `.sif` (`CONTAINER_MODE=singularity`)

One portable file, no registry. Also what you want if you ever need
`--fakeroot` to edit the image.

**On the build server:**

```bash
ssh docker-gpu.wexac.weizmann.ac.il
docker build -f docker/Dockerfile -t phip-pipeline:0.1 docker/
docker images | grep phip-pipeline          # note the IMAGE ID
docker save <IMAGE_ID> -o phip.tar
singularity build --sandbox phip docker-archive://phip.tar
singularity build --fakeroot phip.sif phip
scp phip.sif login1:~/
```

**On Wexac**, the module must be loaded before `singularity` is on PATH:

```bash
bsub -q short -R rusage[mem=16000] 'ml Singularity; singularity exec --bind /home/labs/yourlab:/home/labs/yourlab ~/phip.sif Rscript scripts/04_enrich.R'
```

Set in `config/params.sh`:

```sh
CONTAINER_MODE="singularity"
SIF="$HOME/phip.sif"
BIND_PATHS="/home/labs/yourlab:/home/labs/yourlab"
```

---

## Route C, no build server access

Push the image to Docker Hub from your own machine and let LSF or Singularity
pull it straight from there:

```bash
bsub -q short -app docker-cpu -env LSB_CONTAINER_IMAGE=<youruser>/phip-pipeline:0.1 ...
```

or, converting on the fly (cached under `~/.apptainer/cache/oci-tmp/`):

```bash
bsub -q interactive -Is singularity shell docker://<youruser>/phip-pipeline:0.1
```

---

## Housekeeping

The build server fills up. When you're done testing:

```bash
ssh docker-gpu.wexac.weizmann.ac.il
singularity cache clean
docker image prune
```

## GPU, if you add a structure-prediction stage later

```bash
bsub -q short-gpu -gpu num=1 -app nvidia-gpu -env LSB_CONTAINER_IMAGE=... ./gpu.sh
```

Singularity equivalent needs `--nv`:

```bash
bsub -q short-gpu -gpu num=1 'ml Singularity; singularity exec --nv phip.sif ...'
```

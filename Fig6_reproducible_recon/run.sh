#!/bin/bash
#Copyright 2024. TU Graz. Institute of Biomedical Imaging.
#Author: Moritz Blumenthal

set -eu
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )

export BART_COMPAT_VERSION="v0.9.00"

echo "BART version: $(bart version)"

DDIR=$SCRIPT_DIR/../Data
ODIR=$SCRIPT_DIR/out/
mkdir -p $ODIR

if [ ! -f $SCRIPT_DIR/../Data/ksp_fully.cfl ] ; then

	LINK_TO_ZENODO=???
	wget -o$DDIR/ksp_fully.hdr $LINK_TO_ZENODO/ksp_fully.hdr
	wget -o$DDIR/ksp_fully.cfl $LINK_TO_ZENODO/ksp_fully.cfl
fi

# Working directory for temporary files will be deleted on exit
WORKDIR=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR" || exit

# ## Generate Subsampling Pattern and Undersample k-Space

bart poisson -C24 -Y420 -y4 -Z24 tmp
bart slice 2 0 tmp pat

bart repmat 0 420 pat tmp
cfl2png -x 1 -y 0 tmp $ODIR/pat

bart fmac $DDIR/ksp_fully pat ksp_raw

# ## Prewhiten k-Space Data
# Noise is extracted from a region outside the FOV. The same noise from the undersampled k-space data is used to whiten the the fully sampled and the undersampled k-space data.

bart fft -i -u 1 ksp_raw tmp
bart slice 0 50 tmp tmp
bart transpose 1 0 tmp noise

bart whiten -n ksp_raw noise ksp_white

bart whiten -n $DDIR/ksp_fully noise ksp_fully_white

# ## Fully Sampled Reference Reconstruction

bart ecalib -m1 ksp_fully_white col

bart fft -i -u 3 ksp_fully_white cim_os
bart fmac -C -s8 cim_os col img
bart flip 1 img img

bart resize -c 0 320 1 320 img img
bart rss 8 img $ODIR/rss

cfl2png -A -u800 -x1 -y0 $ODIR/rss $ODIR/rss

# ## Coil Compression to 12 Virtual Coils

bart cc -p 12 ksp_white ksp_cc

# ## Remove Frequency Oversampling

bart fft -u -i 1 ksp_cc tmp
bart resize -c 0 320 tmp tmp
bart fft -u 1 tmp ksp_nos

# ## Reference Reconstruction of Undersampled k-Space Data

POST=_rec

# Estimate Coil Sensitivities
bart ecalib -m1 ksp_cc tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i300 ksp_nos col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST

# ## Missing Information on Number of Iterations

POST=_iter

# Estimate Coil Sensitivities
bart ecalib -m1 ksp_cc tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i100 ksp_nos col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss_rec $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST

# ## Missing Information on Optimization Algorithm (not shown)

POST=_algo

# Estimate Coil Sensitivities
bart ecalib -m1 ksp_cc tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i300 -m ksp_nos col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss_rec $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST

# ## Missing Information on Number of Virtual Coils

POST=_cc

bart cc -p 8 ksp_white ksp_cc2

# Remove frequency oversampling
bart fft -u -i 1 ksp_cc2 tmp
bart resize -c 0 320 tmp tmp
bart fft -u 1 tmp ksp_nos2

# Estimate Coil Sensitivities
bart ecalib -m1 ksp_cc2 tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i300 ksp_nos2 col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss_rec $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST

# ## Missing Information on ESPIRiT Parameter (Not shown)

POST=_coils

# Estimate Coil Sensitivities
# WARNING: on Debian 12 the standard OpenBLAS library triggers a segfault in the SVD
#          (c.f. https://github.com/OpenMathLib/OpenBLAS/issues/5000)
#          As a workaround, the BLAS backend of BART can be changed with
#          $ sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu
#          and selecting the slower reference BLAS implementation.

bart ecalib -a -m1 ksp_cc tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i300 ksp_nos col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss_rec $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST



# ## Missing Information on Normalization of k-Space Data

POST=_normalize

# Estimate Coil Sensitivities
bart ecalib -m1 ksp_cc tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i300 -w1. ksp_nos col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss_rec $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST



# ## Missing Information on Prewhitening

POST=_white

bart cc -p 12 ksp_raw ksp_cc2

# Remove frequency oversampling
bart fft -u -i 1 ksp_cc2 tmp
bart resize -c 0 320 tmp tmp
bart fft -u 1 tmp ksp_nos2

# Estimate Coil Sensitivities
bart ecalib -m1 ksp_cc2 tmp
bart resize -c 0 320 tmp col

# Perform PICS reconstruction
bart pics -S -RW:3:0:0.001 -i300 ksp_nos2 col img
bart flip 1 img img

# Take magnitude
bart rss 0 img rss$POST

# Remove phase oversampling
bart resize -c 0 320 1 320 rss$POST rss$POST

bart nrmse -s rss$POST $ODIR/rss_rec
bart scale 960703.75 rss$POST $ODIR/rss$POST

# Compute difference
bart saxpy -- -1 $ODIR/rss_rec $ODIR/rss$POST $ODIR/diff$POST
bart rss 0 $ODIR/diff$POST $ODIR/diff$POST

cfl2png -A -u800 -x1 -y0 $ODIR/rss$POST $ODIR/rss$POST
cfl2png -A -u40  -x1 -y0 $ODIR/diff$POST $ODIR/diff$POST

# # Generate Final Figure

cd $ODIR
cp ../Fig_6_template.svg Fig_6_reproducible_recon.svg
inkscape --export-filename=../Fig_6_reproducible_recon.pdf Fig_6_reproducible_recon.svg

# %%




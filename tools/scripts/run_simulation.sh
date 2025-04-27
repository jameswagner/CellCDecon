#!/bin/bash
# Script to generate synthetic methylation data and run CellCDecon on it

# Set parameters
SAMPLES=50
PROBES=200
CELL_TYPES=3
OUTPUT_DIR="data/synthetic_data"
OUTPUT_FILE="$OUTPUT_DIR/synthetic_methylation.txt"
NOISE=0.03
CONCENTRATION=0.5  # Lower values create more extreme cell type distributions

# Define paths relative to the project root
PROJECT_ROOT="$(dirname "$(dirname "$(dirname "$(readlink -f "$0")")")")"
GENERATE_SCRIPT="$PROJECT_ROOT/tools/data_generation/generate_test_data.py"
EVALUATE_SCRIPT="$PROJECT_ROOT/tools/analysis/evaluate_results.py"
CELLCDECON_EXE="$PROJECT_ROOT/build/cellcdecon_app"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

echo "Step 1: Generate synthetic methylation data"
python3 $GENERATE_SCRIPT \
  --samples $SAMPLES \
  --probes $PROBES \
  --cell-types $CELL_TYPES \
  --output $OUTPUT_FILE \
  --noise $NOISE \
  --concentration $CONCENTRATION

echo ""
echo "Step 2: Run CellCDecon on the synthetic data"
SEED=$(date +%s)
echo "Running: $CELLCDECON_EXE -k $CELL_TYPES -n $SAMPLES -f $OUTPUT_FILE -c 1 -g 0.005 -s $SEED"
$CELLCDECON_EXE -k $CELL_TYPES -n $SAMPLES -f $OUTPUT_FILE -c 1 -g 0.005 -s $SEED

echo ""
echo "Step 3: Check the results"
echo "CellCDecon output files should be in the $OUTPUT_DIR directory with names:"
echo "- $OUTPUT_FILE.k${CELL_TYPES}.seed${SEED}.gamma0.005.SUMSQ.w       (inferred cell type compositions)"
echo "- $OUTPUT_FILE.k${CELL_TYPES}.seed${SEED}.gamma0.005.SUMSQ.meanvar (inferred methylation values)"
echo "- $OUTPUT_FILE.k${CELL_TYPES}.seed${SEED}.gamma0.005.SUMSQ.resid   (residuals after deconvolution)"
echo ""
echo "True values for comparison:"
echo "- $OUTPUT_FILE.cellcomp (true cell type compositions)"
echo "- $OUTPUT_FILE.means    (true methylation values)"

echo ""
echo "Step 4: Evaluate results"
python3 $EVALUATE_SCRIPT --base-file $OUTPUT_FILE --output-dir "$OUTPUT_DIR/evaluation_results"

echo ""
echo "Evaluation complete. Results saved to $OUTPUT_DIR/evaluation_results" 
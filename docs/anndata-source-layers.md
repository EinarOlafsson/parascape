AnnData source layers
====================

The original-record browser accepts separately verified AnnData HDF5 layers in
CSR, CSC or dense storage. Each descriptor identifies the original matrix path.
Rows remain observations and columns remain variables. Alternative layers and
files retain their own descriptors; they are not independent biological replicates.
The source selector shows the original matrix path (for example `X` or
`layers/spliced`), and the source filter also searches layer names and quantities.

The binary-source receipt must cover every stored value and its order, derived
coordinates, recovery of the original dtype bits, and the original axes. Dense
layers must retain every cell, including explicit zero and nonfinite values when
the producer can prove exact recovery. A failed roundtrip requires a recorded gap,
not coercion or replacement. Sparse absent positions remain absent. Decimal
strings are derived representations; the binary source has no original text
tokens. Decimal NaN strings do not encode the original NaN payload.

The source producer must retain axis metadata with its original types and source
categorical declarations. Integer category codes are not measured abundances.
Neither a layer name such as `matrix`, `absolute`, `spliced` or `velocity`, nor the
filename `raw`, proves its units or suitability for count normalization.

This contract does not establish biological quality, stage, strain, sample
independence, gene equivalence across versions or fair RNA/protein comparability.
Those require source metadata review and separate validation. A layer remaining in
the source inventory is not evidence that it has been processed or released.

Format reference: [AnnData on-disk specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
Older compound-dataframe and `h5sparse_format` files predate the modern encoding
attributes and require explicit legacy handling rather than assumed modern fields.

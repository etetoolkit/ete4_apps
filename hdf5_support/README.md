# Notes

## performance issues

### Alignment

  - Sequence stored in numeric format (int8) may be more efficient and allow hack for generating pseudo consensus sequences.
  - chunk size is very relevant parameter 

__TODO__:
- benchmark chunk size **in both dimensions**
- benchmark random access with
  - sequence length from 10K to 10M
  - tree length from 1K to 1M (1M leaves x 10M bp will not be possible :P)

### Tree depth and complexity:

h5nodes stored using hdf5 group hierarchy, each group is a Tree node (_may be overkill..._).

With `h5py.File(..., libver='latest')` tree size is no longer a problem (tree with 100K leaves created in less than a minute vs **1 hour for 10times less**)

__TODO__: understand why alignment size matters in terms of performance when creating h5Tree

__Other solutions (kind of deprecated as its so fast with latest libver)__:
 - implement tree by level (distance from the root in terms of number of intermediate nodes)
    - GOOD: all nodes in a given level will be available at the same time for fast transversal query
    - BAD: adds a layer of complexity for loading data
    - BAD: harder to reorganize trees

### Potentially slow or difficult reordering

__TODO__: 
  - implement tree rearrangements
  - benchmark on the fly rearrangements versus flush at the end of session (like a `Tree.save_h5Tree`) 

__Solutions__:
 - only change ete.Tree object and `flush` these changes to hdf5 once done, or just for saving
 - deal with it... perhaps only suitable for fixed topology.
 - forget about the tree, only use hdf5 to store other kind of data like alignment etc....

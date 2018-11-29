# tidytree 0.2.0

+ compatible with `tibble` v = 2.0.0 (2018-11-29, Thu)
  - change `as_data_frame` method to `as_tibble` since `as_data_frame` was deprecated in `tibble` and not exported as generics
  
# tidytree 0.1.9

+ `as_data_frame.phylo` works with `phylo$root.edge` (2018-06-13, Wed)

# tidytree 0.1.8

+ force `get_tree_data(treedata)$node` to be integer (2018-04-19, Thu)

# tidytree 0.1.7

+ `get.data`, `[` and `[[` methods (2018-02-26, Mon)

context("tree_subset")

test_that("tree_subset validates selected nodes", {
    tree <- ape::rtree(6)

    expect_error(
        tree_subset(tree, node = "not-a-node"),
        "was not found in the tree",
        fixed = TRUE
    )
})

test_that("tree_subset.treedata passes root_edge to internal subsetter", {
    nwk <- "(((t1:1,t2:1):2,t3:3):4,(t4:1,t5:1):2);"
    tree <- ape::read.tree(text = nwk)
    td <- treedata(phylo = tree)

    with_root_edge <- tree_subset(td, node = "t1", levels_back = 1,
                                  root_edge = TRUE, group_node = FALSE)
    without_root_edge <- tree_subset(td, node = "t1", levels_back = 1,
                                     root_edge = FALSE, group_node = FALSE)

    expect_equal(with_root_edge@phylo$root.edge, 6)
    expect_null(without_root_edge@phylo$root.edge)
})

test_that("as.phylo.tbl_tree fails clearly for malformed tbl_tree objects", {
    x <- tibble::tibble(parent = c(1, 1), node = c(2, 3))
    class(x) <- c("tbl_tree", class(x))

    expect_error(
        as.phylo(x),
        "missing column",
        fixed = TRUE
    )
})

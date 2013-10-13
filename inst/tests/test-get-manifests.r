# test code for get_manifests()

A = matrix(1:12, 4, 3)
colnames(A) = c("x1", "x2", "x3")
rownames(A) = 1:4
a_blocks = list(3, 2, 1, 2:3)
get_manifests(A, a_blocks)

B = as.data.frame(A)
b_blocks = list(2:3, 3:1, 2:1)
B_df = get_manifests(B, b_blocks)
class(B_df)
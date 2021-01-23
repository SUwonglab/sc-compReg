mfbs <- function(TF.name,
                 elem.name,
                 motif.target.path,
                 motif.name,
                 motif.weight,
                 match2,
                 motif.mat.path,
                 ...) {
    UseMethod("mfbs.default")
}

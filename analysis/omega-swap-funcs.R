# Shared by analysis/13-qcs-hs-omega-swap.R and
# analysis/16-qcs-hs-residual-seesaw-experiments.R

# Shift omega (one delta component, at the mesh nodes) so each region's
# area-weighted mean on the grid moves exactly to a target. Nodes take the
# region of their nearest grid cell and each region's nodes get one shift; the
# projection mixes the two shifts across the HS / QCS boundary, so they're
# solved from the 2 x 2 system of how each moves both region means.
#
# This equalizes mean log density (omega), not mean density: with different
# within-region variances, a density contrast can remain after flattening
# (rep_truth()'s contrast_density)
make_region_shifter <- function(mesh, grid) {
  xy <- as.matrix(grid[, c("X", "Y")])
  A <- fmesher::fm_basis(mesh$mesh, loc = xy)
  node_hs <- grid$hs[RANN::nn2(xy, mesh$mesh$loc[, 1:2], k = 1)$nn.idx[, 1]]
  region_area <- c(sum(grid$area[grid$hs == 0L]), sum(grid$area[grid$hs == 1L]))
  # rows: QCS, HS
  W <- rbind(
    (grid$hs == 0L) * grid$area / region_area[1],
    (grid$hs == 1L) * grid$area / region_area[2]
  )
  WA <- as.matrix(W %*% A)
  means <- \(omega) as.vector(WA %*% omega)
  # nodes x regions (QCS, HS) indicator; region means move by M %*% shifts
  D <- cbind(node_hs == 0L, node_hs == 1L) * 1
  M <- WA %*% D
  list(
    project = \(omega) as.vector(A %*% omega),
    means = means,
    overall = \(omega) sum(means(omega) * region_area) / sum(region_area),
    shift = \(omega, target) omega + as.vector(D %*% solve(M, target - means(omega))),
    # Shift omega (the 2nd delta component; `base` is the other component's
    # projection on the grid) so each region's log area-weighted mean of
    # exp(base + omega), the log mean density from omega, moves to the
    # coast-wide value: no contrast in biomass per area. Newton steps with the
    # exact Jacobian (nonlinear because of exp and the boundary mixing)
    shift_density = \(omega, base = 0) {
      for (i in seq_len(50L)) {
        e <- exp(base + as.vector(A %*% omega))
        region_mean <- as.vector(W %*% e)
        lme <- log(region_mean)
        if (abs(diff(lme)) < 1e-8) return(omega)
        target <- log(sum(region_mean * region_area) / sum(region_area))
        J <- as.matrix((W * rep(e, each = 2L)) %*% A %*% D) / region_mean
        omega <- omega + as.vector(D %*% solve(J, target - lme))
      }
      stop("shift_density() did not converge")
    }
  )
}

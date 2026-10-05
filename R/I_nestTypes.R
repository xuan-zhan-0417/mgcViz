##########
# Type of each margin of a gamFactory nested effect ("si", "exp", "mgks", "si_nexp" or "plain"),
# read from the name of its bundle (gamFactory stores si$margin[[k]]$bundle_nam, e.g. "bundle_exp").
.nestTypes <- function(si) {
  sub("bundle_", "", vapply(si$margin, "[[", character(1), "bundle_nam"))
}

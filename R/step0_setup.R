ensure_package_installed <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    install.packages(pkg_name)
  }
}

# install packages
ensure_package_installed("tidyverse")
ensure_package_installed("furrr")
ensure_package_installed("randtoolbox")
ensure_package_installed("prodlim")
ensure_package_installed("matlab")
ensure_package_installed("tmvmixnorm")
ensure_package_installed("reshape2")

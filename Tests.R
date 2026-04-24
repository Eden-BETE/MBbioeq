# BLOCK 1 comp
psi0 <- c(CL = 5, V = 50, ka = 1)

omega_init <- matrix(c(
  0.20, 0.05, 0.02,
  0.05, 0.30, 0.04,
  0.02, 0.04, 0.25
), nrow = 3, byrow = TRUE)


res_block <- run_MBTost(
  data_path = "dataset_1.txt",
  model = "1cpt",
  psi0 = psi0,
  omega_init = omega_init,
  error_model = "combined",
  error_init = c(0.1, 0.2),
)

res_bot <- run_MBbot(
  data_path = "dataset_1.txt",
  model = "1cpt",
  psi0 = psi0,
  omega_init = omega_init,
  error_model = "combined",
  error_init = c(0.1, 0.2)
)



#1 compt avec tlag

psi0 <- c(CL = 4.5, V = 52, ka = 1.1, tlag = 0.5)

omega_init <- c(0.15, 0.25, 0.30, 0.10)   # 4 variances = diagonal

error_init <- c(0.12, 0.25)

res_1cpt_tlag_diag <- run_MBTost(
  data_path = "dataset_1.txt",
  model = "1cpt_tlag",
  psi0 = psi0,
  omega_init = omega_init,
  error_model = "combined",
  error_init = error_init,
)

resbot_1cpt_tlag_diag <- run_MBbot(
  data_path = "dataset_1.txt",
  model = "1cpt_tlag",
  psi0 = psi0,
  omega_init = omega_init,
  error_model = "combined",
  error_init = error_init,
)

# 1 compt avec tlag avec une valeur non corrélée
psi0 <- c(CL = 5, V = 50, ka = 1, tlag = 0.5)
omega_init <- matrix(c(
  0.20, 0,    0,    0,      # CL
  0,    0.30, 0,    0,      # V
  0,    0,    0.25, 0.05,   # ka
  0,    0,    0.05, 0.20    # tlag
), nrow = 4, byrow = TRUE)
res_block_ka_tlag <- run_tost(
  data_path = "dataset_1.txt",
  model = "1cpt_tlag",
  psi0 = psi0,
  omega_init = omega_init,
  covariance_structure = "full",   # car on fournit une matrice complète
  error_model = "combined",
  error_init = c(0.1, 0.2),
)


# Test des plots
plot_spaghetti("dataset_1.txt")
plot_tost(res_block)
plot_bot(res_bot)

import numpy as np

# FIXO: coordenadas da BEPA
X_BEPA = 4229786.5324
Y_BEPA = -4771063.6244
Z_BEPA = -161510.2200

def coord_bepa(comp):
    if comp == "X":
        return X_BEPA
    if comp == "Y":
        return Y_BEPA
    return Z_BEPA

# Dados observados (baseline)
# (i, j, dX, dY, dZ, sigma_mm)
# d = coord(j) - coord(i)
baselines = [
    ("BEPA","M01",  7849.9087,  3085.7272,  1505.4325, 13.80),
    ("M01","M02",   5118.6087,   576.9179,  3131.5131, 16.68),
    ("M02","M03",  -6554.1662,  4284.0852,   223.2862, 13.94),
    ("BEPA","M02", 12968.5476,  3662.5475,  4636.9275, 17.78),
    ("M03","BEPA", -6414.3606, -7946.6716, -4860.2321, 20.42),
]

# Ordem das incógnitas: M01, M02, M03
idx = {
    "M01_X":0, "M01_Y":1, "M01_Z":2,
    "M02_X":3, "M02_Y":4, "M02_Z":5,
    "M03_X":6, "M03_Y":7, "M03_Z":8,
}

# Montagem de A, L, sigma, labels
# (A igual ao seu desenho)
A_rows, L_list, sig_list, labels = [], [], [], []

def add_eq(i, j, d, sigma_m, comp):
    row = np.zeros(9, float)
    label = f"{i}-{j} Δ{comp}"

    # Caso 1: BEPA -> Mk
    # d = Mk - BEPA  => Mk = d + BEPA
    if i == "BEPA" and j != "BEPA":
        row[idx[f"{j}_{comp}"]] = 1.0
        L = d + coord_bepa(comp)

    # Caso 2: Mk -> BEPA (ESTE é o caso do seu desenho)
    # d = BEPA - Mk  => -Mk = d - BEPA
    elif j == "BEPA" and i != "BEPA":
        row[idx[f"{i}_{comp}"]] = -1.0   # <- aqui fica -1, igual ao seu A
        L = d - coord_bepa(comp)         # <- L coerente com isso

    # Caso 3: Mi -> Mj
    # d = Mj - Mi => -Mi + Mj = d
    else:
        row[idx[f"{i}_{comp}"]] = -1.0
        row[idx[f"{j}_{comp}"]] =  1.0
        L = d

    A_rows.append(row)
    L_list.append(L)
    sig_list.append(sigma_m)
    labels.append(label)

for (i,j,dX,dY,dZ,sigma_mm) in baselines:
    s = sigma_mm/1000.0
    add_eq(i,j,dX,s,"X")
    add_eq(i,j,dY,s,"Y")
    add_eq(i,j,dZ,s,"Z")

A = np.vstack(A_rows)
L = np.array(L_list, float)
sigma = np.array(sig_list, float)

# Ajustamento + IDS (1 outlier por vez)
critico = 3.29
ativos = np.arange(len(L))
removidos = []
it = 1

def ajusta(A_sub, L_sub, sigma_sub):
    # 1) Σ_L: matriz de covariância das observações (descorrelacionadas)
    Sigma_L = np.diag(sigma_sub**2)               # em m²

    # 2) P = Σ_L⁻¹
    P = np.diag(1.0/(sigma_sub**2))

    # 3) Ajustamento por mínimos quadrados: x̂ = (Aᵀ P A)⁻¹ Aᵀ P L
    N = A_sub.T @ P @ A_sub                       # N = Aᵀ P A
    n = A_sub.T @ P @ L_sub                       # n = Aᵀ P L
    x = np.linalg.solve(N, n)                     # x̂

    # 4) Resíduos: v = A x̂ − L
    v = A_sub @ x - L_sub                         # em m

    # 5) Construindo H:
    #    L̂ = H L
    #    H = A (Aᵀ P A)⁻¹ Aᵀ P = A N⁻¹ Aᵀ P
    Ninv = np.linalg.inv(N)
    H = A_sub @ Ninv @ A_sub.T @ P                # dimensão m×m

    # 6) Lei geral de propagação:
    #    v = (H − I) L  =>  G = (H − I)
    #    Σ_v = G Σ_L Gᵀ
    I = np.eye(len(L_sub))
    G = (H - I)
    Sigma_v = G @ Sigma_L @ G.T                   # em m²

    # 7) σ_vi e w_i = v_i / σ_vi
    sigma_v = np.sqrt(np.clip(np.diag(Sigma_v), 1e-20, None))  # em m
    w = v / sigma_v                                            # adimensional

    return x, v, sigma_v, w

print("IDS (Data Snooping Iterativo)  |  crítico = 3.29")


while True:
    A_sub = A[ativos,:]
    L_sub = L[ativos]
    s_sub = sigma[ativos]

    x_hat, v, sigma_v, w = ajusta(A_sub, L_sub, s_sub)

    abs_w = np.abs(w)
    k_sub = int(np.argmax(abs_w))
    max_w = float(abs_w[k_sub])

    print(f"\nIteração {it}")
    print(f"Maior |w| = {max_w:.4f}")

    if max_w <= critico:
        print("=> Nenhum outlier. Rede consistente.")
        break

    k_global = int(ativos[k_sub])
    print("=> OUTLIER encontrado:")
    print(f"   Obs: {labels[k_global]}")
    print(f"   w  : {w[k_sub]:.4f}")

    removidos.append((k_global, labels[k_global], float(w[k_sub])))
    ativos = np.delete(ativos, k_sub)
    it += 1

# Impressão final
print("RELATÓRIO FINAL")
print("-"*70)

if len(removidos) == 0:
    print("\nNenhuma observação removida.")
else:
    print("\nOutliers removidos (1 por iteração):")
    for k, lab, wv in removidos:
        print(f"- índice {k:02d} | {lab} | w = {wv:.4f}")
np.set_printoptions(suppress=True)

print(f"M01: X={x_hat[idx['M01_X']]:.4f}  Y={x_hat[idx['M01_Y']]:.4f}  Z={x_hat[idx['M01_Z']]:.4f}")
print(f"M02: X={x_hat[idx['M02_X']]:.4f}  Y={x_hat[idx['M02_Y']]:.4f}  Z={x_hat[idx['M02_Z']]:.4f}")
print(f"M03: X={x_hat[idx['M03_X']]:.4f}  Y={x_hat[idx['M03_Y']]:.4f}  Z={x_hat[idx['M03_Z']]:.4f}")

print(f"\nTotal inicial: {len(L)} | Total final: {len(ativos)}")
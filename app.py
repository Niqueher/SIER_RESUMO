import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Modelo SEIRD-M
def seirdm(y, t, beta, sigma, alfa, gamma, gammaq, mu, muq, N):
    S, E, I, Q, R, D = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma * E
    dIdt = sigma * E - gamma * I - mu * I - alfa * I
    dQdt = alfa * I - (gammaq + muq) * Q
    dRdt = gamma * I + gammaq * Q
    dDdt = mu * I + muq * Q
    return [dSdt, dEdt, dIdt, dQdt, dRdt, dDdt]

st.title("Simulador SEIRD-M")

# Entradas
N = st.number_input("População Total", value=200_000)
beta = st.slider("Taxa de Transmissão (β)", 0.1, 2.0, 1.2)
alfa = st.slider("Taxa de Quarentena (α)", 0.0, 1.0, 0.2)
sigma = 1/5.2
gamma = 1/14
gammaq = 1/2
mu = 0.25
muq = 0.125

# Condições iniciais
S0, E0, I0, Q0, R0, D0 = N - 1, 0, 1, 0, 0, 0
y0 = [S0, E0, I0, Q0, R0, D0]
t = np.linspace(0, 180, 180)

# Resolver ODE
sol = odeint(seirdm, y0, t, args=(beta, sigma, alfa, gamma, gammaq, mu, muq, N))
S, E, I, Q, R, D = sol.T

# Gráfico
fig, ax = plt.subplots()
ax.plot(t, S, label="Suscetíveis")
ax.plot(t, E, label="Expostos")
ax.plot(t, I, label="Infectados")
ax.plot(t, Q, label="Quarentena")
ax.plot(t, R, label="Recuperados")
ax.plot(t, D, label="Mortos")
ax.set_xlabel("Dias")
ax.set_ylabel("Número de Pessoas")
ax.set_title("Dinâmica da Doença")
ax.legend()
st.pyplot(fig)

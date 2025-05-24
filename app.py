import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

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
beta = st.text_input("Taxa de Transmissão (β)", "")
alfa = st.text_input("Taxa de Quarentena (α)", "")
sigma = st.text_input("Tempo de Incubação em Dias (σ)", "")
gamma = st.text_input("Tempo de Recuperação em Dias (σ)", "")
gammaq = st.text_input("Tempo de incubação em dias de Indivíduos em Quarentena (σq)", "")
mu = st.text_input("Taxa de Letalidade (μ)", "")
muq = st.text_input("Taxa de Letalidade de Indivíduos em Quarentena (μq)", "")
time = st.text_input("Número máximo de dias", "")

# Condições iniciais
S0, E0, I0, Q0, R0, D0 = N - 1, 0, 1, 0, 0, 0
y0 = [S0, E0, I0, Q0, R0, D0]
t = np.linspace(0, time, time)

# Resolver ODE
sol = odeint(seirdm, y0, t, args=(pd.to_numeric(beta), 1/sigma, alfa, gamma, gammaq, mu, muq, N))
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

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
from matplotlib.animation import FuncAnimation
import io
from PIL import Image
import base64

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
gammaq = st.text_input("Tempo de Recuperação em Dias de Indivíduos em Quarentena (σq)", "")
mu = st.text_input("Taxa de Letalidade (μ)", "")
muq = st.text_input("Taxa de Letalidade de Indivíduos em Quarentena (μq)", "")
days = pd.to_numeric(st.text_input("Tempo Máximo em Dias", ""))
years = pd.to_numeric(st.text_input("Tempo Máximo em Anos", ""))

# Condições iniciais
S0, E0, I0, Q0, R0, D0 = N - 1, 0, 1, 0, 0, 0
y0 = [S0, E0, I0, Q0, R0, D0]

if not np.isnan(days):
    t = np.linspace(0, days, days)
else:
    t = np.linspace(0, years, years)

# Resolver ODE
sol = odeint(seirdm, y0, t, args=(pd.to_numeric(beta), 1/pd.to_numeric(sigma), pd.to_numeric(alfa), 1/pd.to_numeric(gamma), 1/pd.to_numeric(gammaq), 
                                  pd.to_numeric(mu), pd.to_numeric(muq), pd.to_numeric(N)))
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


# Criar figura para animação
fig2, ax2 = plt.subplots()
lines = {
    'S': ax2.plot([], [], label="Suscetíveis")[0],
    'E': ax2.plot([], [], label="Expostos")[0],
    'I': ax2.plot([], [], label="Infectados")[0],
    'Q': ax2.plot([], [], label="Quarentena")[0],
    'R': ax2.plot([], [], label="Recuperados")[0],
    'D': ax2.plot([], [], label="Mortos")[0],
}
ax2.set_xlim(0, time)
ax2.set_ylim(0, N)
ax2.set_xlabel("Dias")
ax2.set_ylabel("Número de Pessoas")
ax2.set_title("Evolução da Doença ao Longo do Tempo")
ax2.legend()

def update(frame):
    for key, data in zip(['S', 'E', 'I', 'Q', 'R', 'D'], [S, E, I, Q, R, D]):
        lines[key].set_data(t[:frame], data[:frame])
    return lines.values()

ani = FuncAnimation(fig2, update, frames=len(t), interval=50, blit=True)

# Salvar como GIF em memória
gif_buffer = io.BytesIO()
ani.save(gif_buffer, writer='pillow', format='gif')
gif_buffer.seek(0)

# Mostrar GIF no Streamlit
st.image(gif_buffer, caption="Evolução da Epidemia (GIF)", use_column_width=True)

# Criar botão de download
b64 = base64.b64encode(gif_buffer.read()).decode()
href = f'<a href="data:application/octet-stream;base64,{b64}" download="seirdm_animacao.gif">📥 Baixar animação em GIF</a>'
st.markdown(href, unsafe_allow_html=True)

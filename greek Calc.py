import math
from scipy.stats import norm

def black_scholes_greeks(S, K, T, r, sigma):
    # S: spot price
    # K: strike price
    # T: time to maturity
    # r: interest rate
    # sigma: volatility of underlying asset

    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = (math.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))

    delta = norm.cdf(d1)
    gamma = norm.pdf(d1) / (S * sigma * math.sqrt(T))
    theta =- (S * norm.pdf(d1) * sigma / (2 * math.sqrt(T))) - r * K * math.exp(-r * T) * norm.cdf(d2)
    vega = S * norm.pdf(d1) * math.sqrt(T)
    rho = K * T * math.exp(-r * T) * norm.cdf(d2)

    return delta, gamma, theta, vega, rho

S = 51577
K = 51577
T = 1
r = 0.05
sigma = 0.25

delta, gamma, theta, vega, rho = black_scholes_greeks(S, K, T, r, sigma)

print(f"Delta: {delta}, Gamma: {gamma}, Theta: {theta}, Vega: {vega}, Rho: {rho}")
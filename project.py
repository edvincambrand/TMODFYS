import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def get_val(row, i):
    list_values = row.split(";")
    value = list_values[i]
    if "," in value:
        decimals = value.split(",")
        return float(".".join(decimals))
    else:
        return float(value)

def SNR_plots():

    tau = [1,2,3,4,5,6]
    coincidences = [245, 460, 586, 646, 660, 658]
    accidentals = [0.7233363, 1.4312394, 2.1582045, 2.865624, 3.615531, 4.31892]
    SNR = [611.5, 459, 307.4210526, 178.4444444, 145.6666667, 116.5]
    
    ylabels = [r"$R_{T&B}$ (Hz)", "$R_{acc}$ (Hz)", "SNR"]
    titles = ["Coincidence rates", "Accidental coincidences", "Signal-to-noise ratio"]
    colors = ["tab:blue", "tab:orange", "tab:purple"]

    for i, list in enumerate([coincidences, accidentals, SNR]):
        plt.figure()
        plt.plot(tau, list, "o", color=colors[i], ms=6)
        plt.plot(tau, list, "-", color=colors[i], lw=3)
        plt.xlabel(r"Coincidence window $\tau$ (ns)",size=12)
        plt.ylabel(ylabels[i],size=12)
        plt.title(titles[i], size=12)
        plt.grid(True)
        plt.tight_layout()
    
    plt.show()

def coherence_plot():
    stage_positions = []
    coincidences = []


    files = ["9_to_12_step_50_minus46um.csv", "9_to_12_step_50_minus36um.csv",
        "9_to_12_step_50_minus26um.csv", "9_to_12_step_50_minus16um.csv",
        "9_to_12_step_50_minus6um.csv", "9_to_12_step_50_4um.csv",
        "9_to_10_step_50_14um.csv", "9_to_12_step_50_24um.csv",
        "9_to_12_step_50_34um.csv", "9_to_12_step_50_44um.csv",
        "9_to_12_step_50_54um.csv", "9_to_12_step_50_64um.csv"]

    for i, file_name in enumerate(files):
        dx = -46 + 10*i
        with open(file_name) as csvfile:
            for k, row in enumerate(csvfile):
                if k == 0:
                    continue
                stage_positions.append(get_val(row, 0) + dx)
                coincidences.append(get_val(row, 1))
    plt.figure()
    plt.plot(stage_positions, coincidences, "bo", ms=2)
    plt.ylabel(r"Coincidence rates $R_{T&B}$ (Hz)",size=12)
    plt.xlabel(r"Stage position ($\mu$m)",size=12)
    plt.title("Estimating coherence length",size=12)
    plt.tight_layout()
    plt.grid(True)
    plt.show()

def interference():
    coincidences = []
    stage_positions = []
    file_name = "integration_1000_step_25_window_5ns.csv"
    with open(file_name) as csvfile:
        for k, row in enumerate(csvfile):
            if k == 0:
                continue
            stage_positions.append(get_val(row, 0))
            coincidences.append(get_val(row, 1))
    
    plt.figure()
    

    def fit(x, a, b, c, d):
        return a * np.sin(b*x + c) + d
    p0 = [1000, 2*np.pi/0.4, 0, 1400]
    par, parcov = curve_fit(f=fit, xdata=stage_positions, ydata=coincidences, p0=p0)
    x = np.linspace(np.min(stage_positions), np.max(stage_positions), 2000)
    y = par[0] * np.sin(par[1]*x + par[2]) + par[3]
    
    
    plt.plot(stage_positions, coincidences, "ro", ms=4)
    plt.plot(x, y, "-", color="tab:blue", lw=2)


    plt.legend(["Measured", rf"Fitted: $\Delta L = {round(2*np.pi/par[1], 4)} \mu$m"], loc=1)
    plt.xlabel(r"Stage position ($\mu$m)", size=12)
    plt.ylabel(r"Coincidence rate $R_{T&B}$ (Hz)",size=12)
    plt.title("Single photon interference", size=12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def time_delay():
    coincidences = []
    dt = []
    file_name = "coincidence_TB.csv"
    with open(file_name) as csvfile:
        for k, row in enumerate(csvfile):
            if k == 0:
                continue
            dt.append(get_val(row, 0))
            coincidences.append(get_val(row, 2))
    plt.figure()
    plt.plot(dt, coincidences, "r-", lw=2)
    plt.title("Delay adjust", size=12)
    plt.xlabel(r"Time delay $\Delta$t (ns)", size=12)
    plt.ylabel(r"Coincidence rates $R_{T&B}$ (Hz)", size=12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def ryans_plots():
    file_names = ["0_90_22.5_0.csv",
                  "0_90_22.5_22.5.csv",
                  "0_90_22.5_45.csv",
                  "0-90-22.5.csv",
                  "90-90-22.5.csv"]
    titles = [f"Interference. angle: {0}",
              f"Interference. angle: {22.5}",
              f"Interference. angle: {45}",
              f"Interference. 0 / 90",
              f"Interference. 90 / 90"]
    for i, file_name in enumerate(file_names):
        pos = []
        coincidences = []
        with open(file_name) as csv:
            for k, row in enumerate(csv):
                if k == 0:
                    continue
                pos.append(get_val(row, 0))
                coincidences.append(get_val(row, 1))
        plt.figure()
        plt.plot(pos, coincidences, "r-", lw=2)
        plt.plot(pos, coincidences, "ro", ms=4)
        plt.xlabel(r"Stage position ($\mu$m)", size=12)
        plt.ylabel(r"Coincidence rate $R_{T&B}$ (Hz)",size=12)
        plt.title(titles[i], size=12)
        plt.tight_layout()
    plt.show()

def visibility():
    stage_positions = []
    max_rates = []
    min_rates = []

    with open("visibility.csv") as file:
        for k, row in enumerate(file):
            if k == 0:
                continue
            stage_positions.append(get_val(row, 0))
            max_rates.append(get_val(row, 2))
            min_rates.append(get_val(row, 3))
    max_rates = np.array(max_rates)
    min_rates = np.array(min_rates)
    visibility = (max_rates - min_rates) / (max_rates + min_rates)

    def sinc2(x, A, k, x0):
        u = k * (x - x0)
        return np.abs( A * np.sinc(u / np.pi) )

    # Initial guess
    p0 = [0.83915639, 0.07153448, -12.09393824]

    param, param_cov = curve_fit(f=sinc2, xdata=stage_positions, ydata=visibility, p0=p0)
    print(param)
    
    plt.figure()
    
    plt.plot(stage_positions, visibility, "o", color="tab:blue", ms=6)


    x = np.linspace(-80,60,5000)
    u = param[1] * (x - param[2])
    y = np.abs(param[0] * np.sinc(u / np.pi))
    plt.plot(x, y, "-", color="tab:orange", lw=2)

    y_thresh = np.max(y)/np.e
    points = []
    for i in range(1, len(y)):
        if (y[i] - y_thresh)*(y[i-1] - y_thresh) < 0:
            points.append((x[i], y[i]))
    x1, y1 = points[0]
    x2, y2 = points[1]
    plt.plot([x1, x2], [y1, y2], "-", color="black", lw=2)
    l_c = np.abs(x2-x1)

    plt.grid(True)
    plt.xlabel(r"Stage position ($\mu$m)", size=12)
    plt.ylabel("Visibility",size=12)
    plt.title("Visibility plot",size=12)
    plt.tight_layout()
    plt.legend(["Measured", "Fitted", rf"$l_c={round(l_c,2)} \mu m$"])
    plt.show()

    
if __name__ == "__main__":
    #SNR_plots()
    #coherence_plot()
    #interference()
    #time_delay()
    #ryans_plots()
    visibility()
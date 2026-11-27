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
        #plt.plot(tau, list, "o", color=colors[i], ms=6)
        #plt.plot(tau, list, "-", color=colors[i], lw=2)
        err = np.sqrt(np.array(list))
        plt.errorbar(tau, list, yerr=err, fmt="o", color=colors[i], capsize=6, elinewidth=2, ms=6)
        plt.xlabel(r"Coincidence window $\tau$ (ns)",size=14)
        plt.ylabel(ylabels[i],size=14)
        plt.title(titles[i], size=14)
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
    stage_positions = np.array(stage_positions) - 14
    err = np.sqrt(np.array(coincidences))
    plt.figure()
    plt.errorbar(stage_positions, coincidences, yerr=err, fmt="o", color="tab:blue", ms=2, capsize=1)
    plt.ylabel(r"Coincidence rates $R_{T&B}$ (Hz)",size=12)
    plt.xlabel(r"Stage position ($\mu$m)",size=12)
    plt.title("Visualizing coherence length",size=12)
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
    
    err = np.sqrt(np.array(coincidences))

    def fit(x, a, b, c, d):
        return a * np.sin(b*x + c) + d
    p0 = [1000, 2*np.pi/0.4, 0, 1400]
    par, parcov = curve_fit(f=fit, xdata=stage_positions, ydata=coincidences, p0=p0)
    x = np.linspace(np.min(stage_positions), np.max(stage_positions), 2000)
    a,b,c,d=par
    y = fit(x,a,b,c,d)
    
    plt.errorbar(stage_positions, coincidences, yerr=err, color="r", ms=3, fmt="o", capsize=3)
    plt.plot(x, y, "-", color="tab:blue", lw=2)

    # Chi2
    pos = np.array(stage_positions)
    ndof = len(pos) - len(par)
    diff2 = np.square(np.array(coincidences) - fit(pos,a,b,c,d))
    var = err**2
    chi2 = np.sum(diff2/var) / ndof

    plt.legend([rf"Fitted: $\Delta L = {round(2*np.pi/par[1], 3)} \mu$m. $\chi^2$/ndof={round(chi2,3)}", "Measured"], loc=1)
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
    err = np.sqrt(np.array(coincidences))
    plt.figure()
    plt.errorbar(dt, coincidences, yerr=err, color="r", fmt="o", ms=5, capsize=3)
    plt.legend([r"Measured. $\sigma = \sqrt{R_{T&B}}$"])
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
        err = np.sqrt(coincidences)
        plt.figure()
        #plt.plot(pos, coincidences, "r-", lw=2)
        plt.errorbar(pos, coincidences, yerr=err, color="red", fmt="o", ms=4, capsize=3)
        plt.grid(True)
        plt.xlabel(r"Stage position ($\mu$m)", size=12)
        plt.ylabel(r"Coincidence rate $R_{T&B}$ (Hz)",size=12)
        plt.title(titles[i], size=12)
        plt.tight_layout()


        p0list=["X", [200, 2*np.pi/0.4, 0, 200], [200, 2*np.pi/0.4, 0, 200], "X", [400, 2*np.pi/0.4, 0, 280]]
        if i in [1,2,4]:
            def sin(x, a, b, c, d):
                return a*np.sin(b*x+c)+d
            p0 = p0list[i]
            par, cov = curve_fit(sin, pos, coincidences, p0=p0, sigma=err, absolute_sigma=True)
            x = np.linspace(np.min(pos), np.max(pos), 2000)
            a,b,c,d = par
            y = sin(x,a,b,c,d)
            diff2 = np.square(np.array(coincidences) - sin(np.array(pos),a,b,c,d))
            ndof = len(pos) - len(par)
            chi2 = np.sum(diff2/err**2) / ndof
            plt.plot(x, y, color="tab:blue", lw=2)
            plt.legend([fr"Fitted. $\chi^2$/ndof={round(chi2,3)}", "Measured."], loc=1)
        else:
            plt.legend(["Measured."])
    plt.show()

def visibility():
    plt.figure()
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
    stage_positions = np.array(stage_positions) + 14
    max_rates = np.array(max_rates)
    min_rates = np.array(min_rates)
    visibility = (max_rates - min_rates) / (max_rates + min_rates)

    dV = np.sqrt(4*(min_rates**2)*max_rates + 4*(max_rates**2)*min_rates) / (max_rates + min_rates)**2
    plt.errorbar(stage_positions, visibility, yerr=dV, fmt="o", ms=4, capsize=3, color="red")


    # Fit function
    def sinc(x, A, k, x0):
        u = k * (x - x0)
        return np.abs( A * np.sinc(u / np.pi) )

    def sinc2(x, A, k, x0):
        u = k * (x - x0)
        return A * np.sinc(u / np.pi) ** 2

    # Initial guess
    p0 = [0.83915639, 0.07153448, -12.09393824]

    par_sinc, _cov = curve_fit(f=sinc, xdata=stage_positions, ydata=visibility, p0=p0)
    par_sinc2, _cov2 = curve_fit(f=sinc2, xdata=stage_positions, ydata=visibility, p0=p0)
    
    x = np.linspace(-60,70,5000)
    A, k, x0 = par_sinc
    A2, k2, x02 = par_sinc2
    y = sinc(x, A, k, x0)
    y2 = sinc2(x, A2, k2, x02)
    plt.plot(x, y, "-", color="tab:blue", lw=3)
    #plt.plot(x, y2, "--", color="tab:orange", lw=2)

    # Chi2:
    pos = np.array(stage_positions)
    var = dV**2
    diff2 = np.square(visibility - sinc(pos,A,k,x0))
    ndof = len(pos) - len(par_sinc)
    chi2 = np.sum(diff2/var) / ndof

    diff_sinc2 = np.square(visibility - sinc2(pos, A2, k2, x02))
    ndof_sinc2 = len(pos) - len(par_sinc2)
    chi2_sinc2 = np.sum(diff_sinc2/var) / ndof_sinc2

    # Coherence length
    y_thresh = np.max(y)/2
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
    plt.legend([rf"|sinc| $\chi^2$/ndof={round(chi2,3)}", rf"$l_c={round(l_c,2)} \mu m$","Measured"], loc=1)
    plt.show()

    
if __name__ == "__main__":
    #time_delay()
    #SNR_plots()
    coherence_plot()
    #interference()
    #ryans_plots()
    visibility()
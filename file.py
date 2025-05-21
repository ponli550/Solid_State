import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure
theta, intensity = np.loadtxt("ITO.ASC", unpack=True)
peaks, _ = find_peaks(intensity, height=np.max(intensity)*0.1, distance=100)
plt.figure(figsize=(12, 5))
plt.plot(theta, intensity, label="Experimental", color="black")
plt.plot(theta[peaks], intensity[peaks], "rx", markersize=10, label="Detected Peaks")
plt.xlabel("2θ (degrees)")
plt.ylabel("Intensity (a.u.)")
plt.legend()
plt.title("Experimental XRD Pattern")
plt.grid(True)
plt.show()
lambda_xrd = 1.5406  
def bragg_d(two_theta):
    return lambda_xrd / (2 * np.sin(np.radians(two_theta / 2)))
experimental_d = bragg_d(theta[peaks])
formatted_theta = [f"{x:.2f}" for x in theta[peaks]]
formatted_d = [f"{x:.3f}" for x in experimental_d]
print("\nExperimental Peaks Analysis:")
print(f"2θ Positions: {formatted_theta}")
print(f"Calculated d-spacings: {formatted_d} Å")
with MPRester("3ZRh3295Q8T1ygRw4JX8HUv4PqD8FTVY") as mpr:
    structure = mpr.get_structure_by_material_id("mp-755027")
xrd_calculator = XRDCalculator(wavelength="CuKa")
theoretical_pattern = xrd_calculator.get_pattern(structure, two_theta_range=(0, 90))
plt.figure(figsize=(12, 5))
plt.plot(theta, intensity, label="Experimental", color="black")
plt.stem(theoretical_pattern.x, theoretical_pattern.y / theoretical_pattern.y.max(), 
         linefmt='r-', markerfmt='ro', basefmt=' ', label="Theoretical Si")
plt.xlabel("2θ (degrees)")
plt.ylabel("Normalized Intensity")
plt.title("Experimental vs Theoretical Silicon XRD")
plt.legend()
plt.grid(True)
plt.show()
def find_nearest(array, value):
    """Find array index closest to a given value"""
    return (np.abs(array - value)).argmin()
print("\nPeak Matching Analysis:")
print("Exp. 2θ\tExp. d\t\tTheo. 2θ\tTheo. d\t\tMatch?")
for i, exp_two_theta in enumerate(theta[peaks]):
    theo_idx = find_nearest(theoretical_pattern.x, exp_two_theta)
    exp_d = experimental_d[i]
    theo_d = theoretical_pattern.d_hkls[theo_idx]
    delta = abs(exp_d - theo_d)

    print(f"{exp_two_theta:.4f}\t\t{exp_d:.4f}\t\t{theoretical_pattern.x[theo_idx]:.4f}\t\t{theo_d:.4f}\t\t{'✓' if delta < 0.6 else '✗'}")
    
peaks, _ = find_peaks(intensity, height=np.max(intensity)*0.1, distance=100)
def find_nearest(array, value):
    return (np.abs(array - value)).argmin()

matched_hkls = []
for exp_theta in theta[peaks]:
    idx = find_nearest(theoretical_pattern.x, exp_theta)
    hkl = theoretical_pattern.hkls[idx][0]["hkl"]  # Extract Miller indices
    hkl_str = "".join(map(str, hkl))  # Convert (1,1,1) → "111"
    matched_hkls.append(f"({hkl_str})")
plt.figure(figsize=(12, 5))
plt.plot(theta, intensity, label="Experimental", color="black")
plt.plot(theta[peaks], intensity[peaks], "rx", markersize=10, label="Detected Peaks")
for i, (x, y) in enumerate(zip(theta[peaks], intensity[peaks])):
    plt.annotate(
        matched_hkls[i],
        xy=(x, y),
        xytext=(0, 10),
        textcoords='offset points',
        fontsize=9,
        color='blue',
        rotation=0  # Optional: rotate for better spacing
    )

plt.xlabel("2θ (degrees)")
plt.ylabel("Intensity (a.u.)")
plt.title("Experimental XRD with Miller Indices")
plt.legend()
plt.grid(True)
plt.show()

import math

exp_data = [
    {"2theta": 21.5390, "d_exp": 4.1224},
    {"2theta": 24.2170, "d_exp": 3.6722},
    {"2theta": 27.9090, "d_exp": 3.1943},
    {"2theta": 30.6130, "d_exp": 2.9180},
    {"2theta": 35.4750, "d_exp": 2.5284},
    {"2theta": 38.0750, "d_exp": 2.3615},
    {"2theta": 41.6370, "d_exp": 2.1674},
    {"2theta": 45.5370, "d_exp": 1.9904},
    {"2theta": 50.9450, "d_exp": 1.7911},
    {"2theta": 60.5650, "d_exp": 1.5276}
]

hkl_list = [
    (2, 1, 1),  
    (2, 1, 1),   
    (2, 2, 1),  
    (2, 2, 2), 
    (2, 2, 0),  
    (3, 1, 1),   
    (3, 3, 2),   
    (3, 3, 1),   
    (4, 3, 2),  
    (4, 4, 4)   
]

lambda_xrd = 1.5406  

sin2theta_over_s_list = []

for i in range(len(exp_data)):
    two_theta = exp_data[i]["2theta"]
    h, k, l = hkl_list[i]
    
    theta_rad = math.radians(two_theta / 2)
    
    sin_theta = math.sin(theta_rad)
    sin2theta = sin_theta ** 2
    
    s = h**2 + k**2 + l**2
    
    sin2theta_over_s = sin2theta / s
    sin2theta_over_s_list.append(sin2theta_over_s)

avg_sin2theta_over_s = sum(sin2theta_over_s_list) / len(sin2theta_over_s_list)

a = lambda_xrd / (2 * math.sqrt(avg_sin2theta_over_s))

print("Lattice Parameter Calculation:")
print("-" * 60)
print(f"Average sin²θ/s = {avg_sin2theta_over_s:.6f}")
print(f"Lattice parameter a = {a:.4f} Å")

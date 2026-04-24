import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.stats import norm, chi2

def box_muller_transform(mu, sigma, n):
    """Генерация выборки методом Бокса-Мюллера"""
    u1 = np.random.rand(n)
    u2 = np.random.rand(n)
    # Стандартное нормальное распределение
    z0 = np.sqrt(-2 * np.log(u1)) * np.cos(2 * np.pi * u2)
    # Масштабирование под заданные mu и sigma
    return z0 * sigma + mu

def calculate_chi_squared(data, mu, sigma, bins_count=10):
    """Расчет критерия Хи-квадрат"""
    counts, bin_edges = np.histogram(data, bins=bins_count)
    n = len(data)
    observed = counts
    expected = []
    
    for i in range(len(bin_edges) - 1):
        # Вероятность попадания в интервал по теоретической функции
        p = norm.cdf(bin_edges[i+1], mu, sigma) - norm.cdf(bin_edges[i], mu, sigma)
        expected.append(p * n)
    
    expected = np.array(expected)
    # Избегаем деления на ноль, если ожидаемое значение слишком мало
    expected = np.where(expected < 0.1, 0.1, expected)
    
    chi_val = np.sum((observed - expected)**2 / expected)
    # Степени свободы: k - 1 - p (где p=2, так как оцениваем среднее и дисперсию)
    df = bins_count - 1 - 2
    if df <= 0: df = 1
    critical_val = chi2.ppf(0.95, df)
    
    return chi_val, critical_val

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Моделирование: Бокс-Мюллер и Хи-квадрат")
        
        # Панель управления (Ввод данных)
        self.input_frame = tk.Frame(root)
        self.input_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)
        
        tk.Label(self.input_frame, text="Мат. ожидание (μ):").grid(row=0, column=0)
        self.entry_mu = tk.Entry(self.input_frame)
        self.entry_mu.insert(0, "0")
        self.entry_mu.grid(row=0, column=1)
        
        tk.Label(self.input_frame, text="Дисперсия (σ²):").grid(row=0, column=2)
        self.entry_var = tk.Entry(self.input_frame)
        self.entry_var.insert(0, "1")
        self.entry_var.grid(row=0, column=3)
        
        tk.Label(self.input_frame, text="Кол-во элементов (N):").grid(row=0, column=4)
        self.entry_n = tk.Entry(self.input_frame)
        self.entry_n.insert(0, "1000")
        self.entry_n.grid(row=0, column=5)
        
        self.btn_run = tk.Button(self.input_frame, text="Рассчитать", command=self.update_plot, bg="lightblue")
        self.btn_run.grid(row=0, column=6, padx=10)

        # Область для графика
        self.fig, self.ax = plt.subplots(figsize=(8, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Панель результатов (Текст под графиком)
        self.res_frame = tk.Frame(root, bg="white", bd=2, relief=tk.SUNKEN)
        self.res_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=10)
        
        self.lbl_mean = tk.Label(self.res_frame, text="Среднее: -", font=("Arial", 11), bg="white")
        self.lbl_mean.pack(anchor="w")
        
        self.lbl_var = tk.Label(self.res_frame, text="Дисперсия: -", font=("Arial", 11), bg="white")
        self.lbl_var.pack(anchor="w")
        
        self.lbl_chi = tk.Label(self.res_frame, text="Хи-квадрат: -", font=("Arial", 12, "bold"), bg="white")
        self.lbl_chi.pack(anchor="w")

    def update_plot(self):
        try:
            mu_input = float(self.entry_mu.get())
            var_input = float(self.entry_var.get())
            n_input = int(self.entry_n.get())
            sigma_input = np.sqrt(var_input)
            
            if var_input <= 0 or n_input <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Ошибка", "Введите корректные числа (Дисперсия и N должны быть > 0)")
            return

        # Генерация выборки
        sample = box_muller_transform(mu_input, sigma_input, n_input)
        
        # Расчет эмпирических данных
        emp_mean = np.mean(sample)
        emp_var = np.var(sample)
        
        err_mean = abs(emp_mean - mu_input) / abs(mu_input) * 100 if mu_input != 0 else abs(emp_mean) * 100
        err_var = abs(emp_var - var_input) / var_input * 100

        # Отрисовка
        self.ax.clear()
        # Гистограмма эмпирических значений
        count, bins, _ = self.ax.hist(sample, bins=15, density=True, alpha=0.4, color='cornflowerblue', edgecolor='red', label='Эмпирические')
        
        # График реальной плотности (нормальное распределение)
        x = np.linspace(min(sample), max(sample), 100)
        self.ax.plot(x, norm.pdf(x, mu_input, sigma_input), color='green', lw=3, label='Теоретические')
        self.ax.set_title("Сравнение распределений")
        self.ax.legend()
        self.canvas.draw()

        # Расчет Хи-квадрат
        chi_val, crit_val = calculate_chi_squared(sample, mu_input, sigma_input)
        hypothesis_accepted = chi_val <= crit_val
        result_color = "green" if hypothesis_accepted else "red"
        hypo_text = "Принимается" if hypothesis_accepted else "Отвергается"

        # Обновление меток
        self.lbl_mean.config(text=f"Average: {emp_mean:.3f} (error = {err_mean:.1f}%)")
        self.lbl_var.config(text=f"Variance: {emp_var:.3f} (error = {err_var:.1f}%)")
        
        comp_sign = ">" if chi_val > crit_val else "≤"
        self.lbl_chi.config(
            text=f"Chi-squared: {chi_val:.2f} {comp_sign} {crit_val:.2f}  |  Гипотеза: {hypo_text}",
            fg=result_color
        )

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()

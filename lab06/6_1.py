import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class ProbabilityGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Probability Distribution Analysis")
        self.root.geometry("700x600")
        self.root.resizable(False, False)
        
        # Теоретические вероятности
        self.probs = [0.2, 0.2, 0.2, 0.1, 0.3]
        self.values = [1, 2, 3, 4, 5]
        
        # GUI Elements
        self.create_widgets()
        
    def create_widgets(self):
        # Левая панель с вероятностями
        left_frame = tk.Frame(self.root)
        left_frame.pack(side=tk.LEFT, padx=20, pady=20)
        
        tk.Label(left_frame, text="Prob 1", font=("Arial", 10)).grid(row=0, column=0, pady=2)
        tk.Label(left_frame, text="Prob 2", font=("Arial", 10)).grid(row=1, column=0, pady=2)
        tk.Label(left_frame, text="Prob 3", font=("Arial", 10)).grid(row=2, column=0, pady=2)
        tk.Label(left_frame, text="Prob 4", font=("Arial", 10)).grid(row=3, column=0, pady=2)
        tk.Label(left_frame, text="Prob 5", font=("Arial", 10)).grid(row=4, column=0, pady=2)
        
        # Поля ввода вероятностей
        self.prob_entries = []
        for i in range(5):
            entry = tk.Entry(left_frame, width=10, justify='center')
            entry.insert(0, str(self.probs[i]))
            entry.grid(row=i, column=1, padx=10, pady=2)
            self.prob_entries.append(entry)
        
        # Правая панель
        right_frame = tk.Frame(self.root)
        right_frame.pack(side=tk.RIGHT, padx=20, pady=20, fill=tk.BOTH, expand=True)
        
        # Количество экспериментов
        exp_frame = tk.Frame(right_frame)
        exp_frame.pack(pady=10)
        tk.Label(exp_frame, text="Number of experiments", font=("Arial", 10)).pack(side=tk.LEFT)
        self.n_exp = tk.StringVar(value="100")
        tk.Entry(exp_frame, textvariable=self.n_exp, width=10, justify='center').pack(side=tk.LEFT, padx=10)
        
        # Кнопка Start
        tk.Button(right_frame, text="Start", command=self.run_analysis,
                 bg='#4CAF50', fg='white', font=("Arial", 11, "bold"),
                 width=15, pady=5).pack(pady=10)
        
        # Результаты
        results_frame = tk.LabelFrame(right_frame, text="Results", font=("Arial", 10, "bold"))
        results_frame.pack(pady=10, fill=tk.X)
        
        self.avg_label = tk.Label(results_frame, text="Average: ", font=("Arial", 10))
        self.avg_label.pack(anchor='w', padx=10, pady=2)
        
        self.var_label = tk.Label(results_frame, text="Variance: ", font=("Arial", 10))
        self.var_label.pack(anchor='w', padx=10, pady=2)
        
        self.chi2_label = tk.Label(results_frame, text="Chi-squared: ", font=("Arial", 10))
        self.chi2_label.pack(anchor='w', padx=10, pady=2)
        
        # График
        self.fig, self.ax = plt.subplots(figsize=(5, 3))
        self.canvas = FigureCanvasTkAgg(self.fig, right_frame)
        self.canvas.get_tk_widget().pack(pady=10, fill=tk.BOTH, expand=True)
    
    def validate_probabilities(self, probs):
        """Проверяет корректность вероятностей"""
        # Проверка на отрицательные значения
        for i, p in enumerate(probs):
            if p < 0:
                messagebox.showerror("Ошибка", 
                                    f"Вероятность {i+1} не может быть отрицательной!\n"
                                    f"Текущее значение: {p:.3f}\n"
                                    f"Пожалуйста, введите неотрицательное число.")
                return False
        
        # Проверка суммы
        prob_sum = sum(probs)
        if abs(prob_sum - 1.0) > 1e-6:
            messagebox.showerror("Ошибка", 
                                f"Сумма вероятностей должна быть равна 1!\n"
                                f"Текущая сумма: {prob_sum:.3f}\n"
                                f"Пожалуйста, исправьте значения вероятностей.")
            return False
        
        return True
    
    def run_analysis(self):
        try:
            # Получаем вероятности
            probs = []
            for entry in self.prob_entries:
                val = float(entry.get())
                probs.append(val)
            
            # Валидация вероятностей
            if not self.validate_probabilities(probs):
                return
            
            self.probs = probs
            n = int(self.n_exp.get())
            
            if n <= 0:
                messagebox.showerror("Ошибка", "Количество экспериментов должно быть положительным числом!")
                return
                
        except ValueError:
            messagebox.showerror("Ошибка", "Пожалуйста, введите корректные числовые значения!")
            return
        
        # Преобразуем в numpy массивы
        probs_array = np.array(self.probs)
        values_array = np.array(self.values)
        
        # Генерация выборки
        cumsum = np.cumsum(probs_array)
        counts = np.zeros(5, dtype=int)
        
        for _ in range(n):
            r = np.random.random()
            for i, cum in enumerate(cumsum):
                if r <= cum:
                    counts[i] += 1
                    break
        
        emp_probs = counts / n
        
        emp_mean = np.sum(counts * values_array) / n
        theor_mean = np.sum(probs_array * values_array)
        mean_error = abs(emp_mean - theor_mean) / theor_mean * 100 if theor_mean != 0 else 0
        
        emp_var = np.sum(counts * (values_array - emp_mean)**2) / n
        theor_var = np.sum(probs_array * (values_array - theor_mean)**2)
        var_error = abs(emp_var - theor_var) / theor_var * 100 if theor_var != 0 else 0
        

        expected = n * probs_array
        
        # Проверка условия применимости критерия
        if np.any(expected < 5):
            warning = "Внимание: некоторые ожидаемые частоты < 5\n"
            warning += "Результаты критерия могут быть ненадежны."
            messagebox.showwarning("Предупреждение", warning)
        
        # Расчет статистики χ²
        chi2_stat = np.sum((counts - expected)**2 / expected)
        

        df = len(self.values) - 1
        critical = chi2.ppf(0.95, df)
        p_value = 1 - chi2.cdf(chi2_stat, df)
        verdict = chi2_stat > critical
        
        # Обновление интерфейса
        self.avg_label.config(text=f"Average: {emp_mean:.3f} (error = {mean_error:.0f}%)")
        self.var_label.config(text=f"Variance: {emp_var:.3f} (error = {var_error:.0f}%)")
        
        if verdict:
            self.chi2_label.config(text=f"Chi-squared: {chi2_stat:.2f} > {critical:.3f} (p={p_value:.3f}) - H0 rejected")
        else:
            self.chi2_label.config(text=f"Chi-squared: {chi2_stat:.2f} <= {critical:.3f} (p={p_value:.3f}) - H0 accepted")
        
        # Обновление графика
        self.ax.clear()
        x = np.arange(len(self.values))
        width = 0.35
        
        bars = self.ax.bar(x, emp_probs, width, label='Empirical', color='skyblue', alpha=0.7)
        self.ax.bar(x + width, probs_array, width, label='Theoretical', color='orange', alpha=0.7)
        
        self.ax.set_xticks(x + width/2)
        self.ax.set_xticklabels(self.values)
        self.ax.set_ylabel('Frequency')
        self.ax.set_xlabel('Values')
        self.ax.set_title(f'Distribution (n={n})')
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
        
        # Добавление значений на столбцы
        for i, (bar, prob) in enumerate(zip(bars, emp_probs)):
            height = bar.get_height()
            self.ax.text(bar.get_x() + bar.get_width()/2., height,
                        f'{prob:.3f}', ha='center', va='bottom', fontsize=9)
        
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = ProbabilityGUI(root)
    root.mainloop()

import tkinter as tk
from tkinter import ttk, messagebox
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class WeatherSimulationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Симуляция погоды (Непрерывное время)")
        
        # --- Параметры интерфейса ---
        self.setup_ui()
        
        # --- Математические параметры ---
        self.states = [0, 1, 2]
        self.names = ["Ясно", "Облачно", "Пасмурно"]
        self.colors = ['#ffcc00', '#aaaaaa', '#555555']
        
    def setup_ui(self):
        # Левая панель настроек
        input_frame = ttk.LabelFrame(self.root, text=" Настройки интенсивностей (Q) ")
        input_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        self.entries = {}
        # Задаем только внедиагональные элементы (i != j)
        pairs = [(0,1), (0,2), (1,0), (1,2), (2,0), (2,1)]
        defaults = [0.4, 0.1, 0.3, 0.4, 0.2, 0.3]

        for idx, (i, j) in enumerate(pairs):
            ttk.Label(input_frame, text=f"λ {i+1}→{j+1}:").grid(row=idx, column=0, sticky=tk.W)
            ent = ttk.Entry(input_frame, width=10)
            ent.insert(0, str(defaults[idx]))
            ent.grid(row=idx, column=1, padx=5, pady=2)
            self.entries[f"{i}{j}"] = ent

        ttk.Label(input_frame, text="Общее время (T):").grid(row=6, column=0, sticky=tk.W, pady=10)
        self.entry_T = ttk.Entry(input_frame, width=10)
        self.entry_T.insert(0, "100")
        self.entry_T.grid(row=6, column=1)

        self.btn_start = ttk.Button(input_frame, text="Запустить", command=self.run_simulation)
        self.btn_start.grid(row=7, column=0, columnspan=2, pady=20)

        # Правая панель графиков
        self.fig, self.ax = plt.subplots(figsize=(5, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def get_q_matrix(self):
        """Формирует матрицу Q на основе ввода, считая диагонали автоматически"""
        Q = np.zeros((3, 3))
        for key, entry in self.entries.items():
            i, j = int(key[0]), int(key[1])
            Q[i, j] = float(entry.get())
        
        # Считаем диагонали: q_ii = -sum(q_ij)
        for i in range(3):
            Q[i, i] = -np.sum(Q[i, :])
        return Q

    def get_stationary_dist(self, Q):
        """Решает систему pi * Q = 0"""
        n = Q.shape[0]
        A = Q.T
        A[-1] = np.ones(n)
        b = np.zeros(n)
        b[-1] = 1.0
        try:
            return np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return np.array([1/3, 1/3, 1/3])

    def run_simulation(self):
        try:
            Q = self.get_q_matrix()
            T_max = float(self.entry_T.get())
        except ValueError:
            messagebox.showerror("Ошибка", "Введите корректные числа")
            return

        stationary_theory = self.get_stationary_dist(Q)
        
        # Матрица вероятностей прыжков (диагонали 0)
        jump_matrix = []
        for i in range(3):
            lambda_out = abs(Q[i, i])
            probs = [(Q[i, j] / lambda_out) if i != j else 0 for j in range(3)]
            jump_matrix.append(probs)

        # --- Инициализация симуляции ---
        t = 0
        ans = list(stationary_theory)
        current_state = ans.index(max(ans))
        stat_time = [0.0, 0.0, 0.0]
        history = [] # Для CSV

        self.btn_start.config(state=tk.DISABLED)

        # Внутренний цикл симуляции (анимация)
        def step():
            nonlocal t, current_state
            
            lambda_out = abs(Q[current_state, current_state])
            dt = - math.log(random.random()) / lambda_out
            
            if t + dt > T_max:
                dt = T_max - t
                stat_time[current_state] += dt
                t = T_max
                finish()
                return

            stat_time[current_state] += dt
            t += dt
            
            # Запись шага
            history.append({"Time": round(t, 4), "State": self.names[current_state]})
            
            # Прыжок (используем заранее посчитанные веса)
            current_state = random.choices(self.states, weights=jump_matrix[current_state])[0]

            # Обновление Pie Chart
            self.update_chart(stat_time)
            
            # Скорость анимации (ms)
            self.root.after(10, step)

        def finish():
            self.btn_start.config(state=tk.NORMAL)
            emp_dist = [s / T_max for s in stat_time]
            
            # Сохранение в CSV
            df = pd.DataFrame(history)
            df.to_csv("simulation_results.csv", index=False)
            
            # Вывод сравнения
            report = "Сравнение распределений:\n\n"
            for i in range(3):
                report += f"{self.names[i]}:\n  Теория: {stationary_theory[i]:.3f}\n  Факт: {emp_dist[i]:.3f}\n\n"
            
            messagebox.showinfo("Результаты", report + "Данные сохранены в simulation_results.csv")

        step()

    def update_chart(self, stat_time):
        self.ax.clear()
        # Чтобы избежать ошибки, если сумма еще 0
        if sum(stat_time) > 0:
            self.ax.pie(stat_time, labels=self.names, colors=self.colors, autopct='%1.1f%%', startangle=140)
        self.ax.set_title("Распределение типов погоды в реальном времени")
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = WeatherSimulationGUI(root)
    root.mainloop()

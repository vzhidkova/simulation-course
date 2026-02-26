import numpy as np
from numba import njit
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from time import time
import threading

class HeatSimulationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Моделирование теплопроводности")
        self.root.geometry("900x800")
        
        # Словарь с материалами
        self.materials = {
            "Медь": {"lam": 389.6, "rho": 8960, "c": 385, "name": "Медь"},
            "Алюминий": {"lam": 237, "rho": 2700, "c": 900, "name": "Алюминий"},
            "Сталь": {"lam": 45, "rho": 7800, "c": 460, "name": "Сталь"}
        }
        
        self.create_widgets()
        
    def create_widgets(self):
        # Основной фрейм для параметров
        params_frame = ttk.LabelFrame(self.root, text="Параметры моделирования", padding=10)
        params_frame.pack(fill="x", padx=10, pady=5)
        
        # Материал
        ttk.Label(params_frame, text="Материал:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.material_var = tk.StringVar(value="Медь")
        material_combo = ttk.Combobox(params_frame, textvariable=self.material_var, 
                                       values=list(self.materials.keys()), state="readonly", width=15)
        material_combo.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        material_combo.bind('<<ComboboxSelected>>', self.update_material_params)
        
        # Параметры материала (отображение)
        ttk.Label(params_frame, text="Теплопроводность:").grid(row=0, column=2, sticky="e", padx=5, pady=5)
        self.lam_var = tk.StringVar(value="389.6")
        ttk.Label(params_frame, textvariable=self.lam_var, width=10).grid(row=0, column=3, sticky="w", padx=5, pady=5)
        
        ttk.Label(params_frame, text="Плотность:").grid(row=1, column=2, sticky="e", padx=5, pady=5)
        self.rho_var = tk.StringVar(value="8960")
        ttk.Label(params_frame, textvariable=self.rho_var, width=10).grid(row=1, column=3, sticky="w", padx=5, pady=5)
        
        ttk.Label(params_frame, text="Теплоемкость:").grid(row=2, column=2, sticky="e", padx=5, pady=5)
        self.c_var = tk.StringVar(value="385")
        ttk.Label(params_frame, textvariable=self.c_var, width=10).grid(row=2, column=3, sticky="w", padx=5, pady=5)
        
        # Геометрические параметры
        ttk.Label(params_frame, text="Длина (м):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.length_var = tk.StringVar(value="0.3")
        ttk.Entry(params_frame, textvariable=self.length_var, width=15).grid(row=1, column=1, sticky="w", padx=5, pady=5)
        
        ttk.Label(params_frame, text="Время моделирования (с):").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.time_var = tk.StringVar(value="2.0")
        ttk.Entry(params_frame, textvariable=self.time_var, width=15).grid(row=2, column=1, sticky="w", padx=5, pady=5)
        
        # Температуры
        temp_frame = ttk.LabelFrame(params_frame, text="Температуры (°C)", padding=5)
        temp_frame.grid(row=3, column=0, columnspan=4, sticky="ew", padx=5, pady=5)
        
        ttk.Label(temp_frame, text="Левая граница:").grid(row=0, column=0, padx=5, pady=2)
        self.t_left_var = tk.StringVar(value="100.0")
        ttk.Entry(temp_frame, textvariable=self.t_left_var, width=10).grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(temp_frame, text="Правая граница:").grid(row=0, column=2, padx=5, pady=2)
        self.t_right_var = tk.StringVar(value="20.0")
        ttk.Entry(temp_frame, textvariable=self.t_right_var, width=10).grid(row=0, column=3, padx=5, pady=2)
        
        ttk.Label(temp_frame, text="Начальная:").grid(row=0, column=4, padx=5, pady=2)
        self.t0_var = tk.StringVar(value="20.0")
        ttk.Entry(temp_frame, textvariable=self.t0_var, width=10).grid(row=0, column=5, padx=5, pady=2)
        
        # Кнопка запуска
        self.run_button = ttk.Button(params_frame, text="Запустить моделирование", command=self.run_simulation)
        self.run_button.grid(row=4, column=0, columnspan=4, pady=10)
        
        # Фрейм для результатов
        results_frame = ttk.LabelFrame(self.root, text="Результаты моделирования", padding=10)
        results_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Текстовое поле для вывода таблиц
        self.results_text = tk.Text(results_frame, height=20, width=100, font=("Courier", 10))
        self.results_text.pack(fill="both", expand=True)
        
        # Скроллбар
        scrollbar = ttk.Scrollbar(self.results_text, command=self.results_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.results_text.config(yscrollcommand=scrollbar.set)
        
    def update_material_params(self, event=None):
        material = self.material_var.get()
        params = self.materials[material]
        self.lam_var.set(str(params["lam"]))
        self.rho_var.set(str(params["rho"]))
        self.c_var.set(str(params["c"]))
        
    def validate_inputs(self):
        try:
            float(self.length_var.get())
            float(self.time_var.get())
            float(self.t_left_var.get())
            float(self.t_right_var.get())
            float(self.t0_var.get())
            float(self.lam_var.get())
            float(self.rho_var.get())
            float(self.c_var.get())
            return True
        except ValueError:
            messagebox.showerror("Ошибка", "Пожалуйста, введите корректные числовые значения")
            return False
    
    @staticmethod
    @njit
    def solve_heat(dt, dx, L, t_target, T_left, T_right, T0, lam, rho, c):
        nx = int(L / dx) + 1
        nt = int(t_target / dt) + 1
        
        T = np.full(nx, T0)
        T[0] = T_left
        T[-1] = T_right
        
        ai = lam / (dx ** 2)
        ci = lam / (dx ** 2)
        
        alphas = np.full(nx, 0.0)
        beths = np.full(nx, 0.0)
        beths[0] = T_left
        
        for t in range(1, nt):
            bi = 2 * lam / (dx ** 2) + (rho * c) / dt
            
            for i in range(1, nx - 1):
                fi = -(c * rho / dt) * T[i]
                alphas[i] = ai / (bi - ci * alphas[i - 1])
                beths[i] = (ci * beths[i - 1] - fi) / (bi - ci * alphas[i - 1])
            
            for i in range(nx - 2, -1, -1):
                T[i] = alphas[i] * T[i + 1] + beths[i]
            
            T[0] = T_left
            T[-1] = T_right
        
        center_temp = T[nx // 2]
        return center_temp
    
    def run_simulation(self):
        if not self.validate_inputs():
            return
        
        # Блокируем кнопку на время расчета
        self.run_button.config(state="disabled")
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(tk.END, "Выполняется моделирование...\n")
        self.root.update()
        
        # Запускаем моделирование в отдельном потоке, чтобы GUI не зависал
        thread = threading.Thread(target=self._run_calculation)
        thread.daemon = True
        thread.start()
    
    def _run_calculation(self):
        try:
            # Получаем параметры
            L = float(self.length_var.get())
            t_target = float(self.time_var.get())
            T_left = float(self.t_left_var.get())
            T_right = float(self.t_right_var.get())
            T0 = float(self.t0_var.get())
            lam = float(self.lam_var.get())
            rho = float(self.rho_var.get())
            c = float(self.c_var.get())
            
            time_steps = [0.1, 0.01, 0.001, 0.0001]
            space_steps = [0.1, 0.01, 0.001, 0.0001]
            
            results = np.zeros((len(space_steps), len(time_steps)))
            times = np.zeros((len(space_steps), len(time_steps)))
            
            # Формируем вывод
            output = "\n" + "="*80 + "\n"
            output += f"МАТЕРИАЛ: {self.material_var.get()}\n"
            output += f"Параметры: λ={lam} Вт/(м·К), ρ={rho} кг/м³, c={c} Дж/(кг·К)\n"
            output += f"Длина: {L} м, Время моделирования: {t_target} с\n"
            output += f"Температуры: Tл={T_left}°C, Tп={T_right}°C, T0={T0}°C\n"
            output += "="*80 + "\n\n"
            
            output += "ТАБЛИЦА 1: Температура в центральной точке после {} с модельного времени (°C)\n".format(t_target)
            output += "-" * 65 + "\n"
            output += "Шаг по времени, с | 0.1      | 0.01     | 0.001    | 0.0001\n"
            output += "-" * 65 + "\n"
            
            for i, dx in enumerate(space_steps):
                output += f"Шаг по пространству {dx:<6} м | "
                for j, dt in enumerate(time_steps):
                    start_time = time()
                    Tc = self.solve_heat(dt, dx, L, t_target, T_left, T_right, T0, lam, rho, c)
                    elapsed_time = time() - start_time
                    results[i, j] = Tc
                    times[i, j] = elapsed_time
                    output += f"{Tc:<8.4f} "
                    
                    # Обновляем GUI через каждые 4 расчета
                    if (i * len(time_steps) + j + 1) % 4 == 0:
                        self.results_text.delete(1.0, tk.END)
                        self.results_text.insert(tk.END, output + "\nВыполняется...")
                        self.root.update()
                
                output += "\n"
            
            output += "-" * 65 + "\n\n"
            
            output += "ТАБЛИЦА 2: Реальное время моделирования (с)\n"
            output += "-" * 65 + "\n"
            output += "Шаг по времени, с | 0.1      | 0.01     | 0.001    | 0.0001\n"
            output += "-" * 65 + "\n"
            
            for i, dx in enumerate(space_steps):
                output += f"Шаг по пространству {dx:<6} м | "
                for j, dt in enumerate(time_steps):
                    output += f"{times[i, j]:<8.6f} "
                output += "\n"
            
            output += "-" * 65 + "\n"
            
            # Выводим результат
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, output)
            
        except Exception as e:
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, f"Ошибка при моделировании: {str(e)}")
        finally:
            # Разблокируем кнопку
            self.run_button.config(state="normal")

if __name__ == '__main__':
    root = tk.Tk()
    app = HeatSimulationGUI(root)
    root.mainloop()

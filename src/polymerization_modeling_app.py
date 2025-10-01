import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import signal
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
import os
import sys

# 设置matplotlib后端和字体
matplotlib.use('TkAgg')
plt.rcParams['axes.unicode_minus'] = False

class PolymerizationKineticsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Polymerization Kinetics Simulator")
        self.root.geometry("1200x800")
        
        # 设置信号处理
        signal.signal(signal.SIGINT, self.handle_interrupt)
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # 默认参数
        self.default_params = {
            'MONOMER': 'MMA',
            'INITIATOR': 'Benzoyl peroxide',
            'TEMPERATURE': 50,
            'KP': 587,
            'KT': 1.37e7,
            'KD': 9.16e-7,
            'INITIATOR_EFFICIENCY': 0.6,
            'CHAIN_TRANSFER_CONSTANT': 1.8,
            'MONOMER_MW': 100.12,
            'INITIATOR_CONC0': 0.005,
            'MONOMER_CONC0': 3.0,
            'TIME_HOURS_MAX': 12.0,
            'TIME_STEP': 0.2,
            'TERMINATION_MODE': 1.0
        }
        
        self.params = self.default_params.copy()
        
        # Arrhenius参数
        self.arrhenius_data = {
            50: {'KP': 587, 'KT': 1.37e7, 'KD': 9.16e-7},
            70: {
                'KP_A': 8.7e5, 'KP_Ea': 4.7e3,
                'KT_A': 9.1e8, 'KT_Ea': 2.7e3,
                'KD_A': 1.58e14, 'KD_Ea': 30.0e3
            }
        }
        
        self.setup_ui()
        self.update_plot()
    
    def handle_interrupt(self, signum, frame):
        """处理Ctrl+C中断"""
        print("\nReceived Ctrl+C. Closing the application...")
        self.root.after(0, self.on_closing)

    def on_closing(self):
        """处理窗口关闭"""
        print("Closing application...")
        plt.close('all')
        self.root.destroy()
        sys.exit(0)
    
    def setup_ui(self):
        # 主框架
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 左侧控制面板
        control_frame = ttk.LabelFrame(main_frame, text="Simulation Parameters", padding=10)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # 右侧图表面板
        plot_frame = ttk.LabelFrame(main_frame, text="Results", padding=10)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # 创建控件
        self.create_controls(control_frame)
        
        # 设置图表区域
        self.setup_plot(plot_frame)
        
        # 理论部分
        self.create_theory_section(control_frame)
    
    def create_controls(self, parent):
        # 温度选择
        ttk.Label(parent, text="Temperature (°C):").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.temp_var = tk.StringVar(value=str(self.params['TEMPERATURE']))
        temp_frame = ttk.Frame(parent)
        temp_frame.grid(row=0, column=1, sticky=tk.W, pady=5)
        ttk.Radiobutton(temp_frame, text="50 °C", variable=self.temp_var, 
                       value="50", command=self.on_parameter_change).pack(side=tk.LEFT)
        ttk.Radiobutton(temp_frame, text="70 °C", variable=self.temp_var, 
                       value="70", command=self.on_parameter_change).pack(side=tk.LEFT)
        
        # 单体浓度
        ttk.Label(parent, text="Initial Monomer Conc. [M]₀ (M):").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.monomer_conc_var = tk.StringVar(value=str(self.params['MONOMER_CONC0']))
        monomer_conc_entry = ttk.Entry(parent, textvariable=self.monomer_conc_var, width=10)
        monomer_conc_entry.grid(row=1, column=1, sticky=tk.W, pady=5)
        monomer_conc_entry.bind('<KeyRelease>', self.on_parameter_change)
        
        # 引发剂浓度
        ttk.Label(parent, text="Initial Initiator Conc. [I]₀ (M):").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.initiator_conc_var = tk.StringVar(value=str(self.params['INITIATOR_CONC0']))
        initiator_conc_entry = ttk.Entry(parent, textvariable=self.initiator_conc_var, width=10)
        initiator_conc_entry.grid(row=2, column=1, sticky=tk.W, pady=5)
        initiator_conc_entry.bind('<KeyRelease>', self.on_parameter_change)
        
        # 终止模式滑块
        ttk.Label(parent, text="Termination Mode (z):").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.termination_var = tk.DoubleVar(value=self.params['TERMINATION_MODE'])
        termination_scale = ttk.Scale(parent, variable=self.termination_var,
                                     from_=1.0, to=2.0, orient=tk.HORIZONTAL,
                                     length=150, command=self.update_termination_label)
        termination_scale.grid(row=3, column=1, sticky=tk.W, pady=5)
        
        # 滑块值显示标签
        self.termination_label_var = tk.StringVar(value=f"{self.termination_var.get():.1f}")
        termination_label = ttk.Label(parent, textvariable=self.termination_label_var, width=5)
        termination_label.grid(row=3, column=2, sticky=tk.W, padx=5)
        
        # 按钮框架
        button_frame = ttk.Frame(parent)
        button_frame.grid(row=4, column=0, columnspan=3, pady=10)
        
        # 更新图表按钮
        ttk.Button(button_frame, text="Update Plot", command=self.update_plot).pack(side=tk.LEFT, padx=5)
        
        # 重置参数按钮
        ttk.Button(button_frame, text="Reset Parameters", command=self.reset_parameters).pack(side=tk.LEFT, padx=5)
        
        # 导出图表按钮
        ttk.Button(button_frame, text="Export Plot (PNG)", command=self.export_plot).pack(side=tk.LEFT, padx=5)
        
        # 参数显示区域
        params_frame = ttk.LabelFrame(parent, text="Current Parameters", padding=5)
        params_frame.grid(row=5, column=0, columnspan=3, sticky=tk.W+tk.E, pady=10)
        
        self.params_text = tk.Text(params_frame, height=8, width=40, font=("Courier", 9))
        self.params_text.pack(fill=tk.BOTH, expand=True)
    
    def update_termination_label(self, value):
        """更新终止模式标签并触发参数变化"""
        self.termination_label_var.set(f"{float(value):.1f}")
        self.on_parameter_change()

    def create_theory_section(self, parent):
        theory_frame = ttk.LabelFrame(parent, text="Theory & Equations", padding=5)
        theory_frame.grid(row=6, column=0, columnspan=3, sticky=tk.W+tk.E, pady=10)
        
        theory_text = tk.Text(theory_frame, height=12, width=40, font=("Arial", 9))
        theory_text.pack(fill=tk.BOTH, expand=True)
        
        theory_content = """
Polymerization Kinetics Model

Key Equations:

1. Monomer Conversion:
   [M]/[M]0 = exp( A × (exp(-kd*t/2) - 1) )
   
   where A = (2kp/kd) × sqrt(f*kd*[I]0/kt)

2. Molecular Weight:
   Mn = (Δ[M] / (f*z*Δ[I])) × MW_monomer

Variables:
- [M]0: Initial monomer conc. (M)
- [I]0: Initial initiator conc. (M)
- kp: Propagation rate constant
- kt: Termination rate constant  
- kd: Initiation rate constant
- f: Initiator efficiency
- z: Chain transfer constant (Termination Mode)

Units:
- Concentrations: M (mol/L)
- Rate constants: vary (1/s, 1/(M·s))
- Molecular weight: g/mol
- Time: hours
"""
        theory_text.insert(tk.END, theory_content)
        theory_text.config(state=tk.DISABLED)
    
    def setup_plot(self, parent):
        self.fig, (self.ax1, self.ax2) = plt.subplots(2, 1, figsize=(8, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def on_parameter_change(self, event=None):
        self.update_parameters_display()
    
    def update_parameters_display(self):
        try:
            temp = int(self.temp_var.get())
            monomer_conc = float(self.monomer_conc_var.get())
            initiator_conc = float(self.initiator_conc_var.get())
            termination_mode = self.termination_var.get()
            
            # 基于温度更新速率常数
            if temp == 50:
                kp = self.arrhenius_data[50]['KP']
                kt = self.arrhenius_data[50]['KT'] * termination_mode
                kd = self.arrhenius_data[50]['KD']
            else:
                R = 1.987  # cal/(mol·K)
                T_K = temp + 273.15
                arr_data = self.arrhenius_data[70]
                kp = arr_data['KP_A'] * np.exp(-arr_data['KP_Ea'] / (R * T_K))
                kt = arr_data['KT_A'] * np.exp(-arr_data['KT_Ea'] / (R * T_K)) * termination_mode
                kd = arr_data['KD_A'] * np.exp(-arr_data['KD_Ea'] / (R * T_K))
            
            params_text = f"""Current Parameters:
Temperature: {temp} °C
[M]_0: {monomer_conc} M
[I]_0: {initiator_conc} M
Termination Mode (z): {termination_mode:.1f}

Rate Constants:
kp: {kp:.3e} 1/(M·s)
kt: {kt:.3e} 1/(M·s)  
kd: {kd:.3e} 1/s

Other:
f: {self.params['INITIATOR_EFFICIENCY']}
MW: {self.params['MONOMER_MW']} g/mol"""
            
            self.params_text.config(state=tk.NORMAL)
            self.params_text.delete(1.0, tk.END)
            self.params_text.insert(tk.END, params_text)
            self.params_text.config(state=tk.DISABLED)
            
        except ValueError:
            pass
    
    def reset_parameters(self):
        """重置所有参数为默认值"""
        self.params = self.default_params.copy()
        self.temp_var.set(str(self.params['TEMPERATURE']))
        self.monomer_conc_var.set(str(self.params['MONOMER_CONC0']))
        self.initiator_conc_var.set(str(self.params['INITIATOR_CONC0']))
        self.termination_var.set(self.params['TERMINATION_MODE'])
        self.update_termination_label(str(self.params['TERMINATION_MODE']))
        self.update_plot()
        messagebox.showinfo("Reset", "Parameters have been reset to default values.")

    def export_plot(self):
        """将当前图表导出为PNG文件"""
        try:
            # 生成基于参数的默认文件名
            temp = int(self.temp_var.get())
            monomer_conc = float(self.monomer_conc_var.get())
            initiator_conc = float(self.initiator_conc_var.get())
            termination_mode = self.termination_var.get()
            
            default_filename = f"polymerization_T{temp}C_M{monomer_conc}_I{initiator_conc}_z{termination_mode:.1f}.png"
            
            # 询问用户保存位置和文件名
            file_path = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
                initialfile=default_filename,
                title="Export Plot As"
            )
            
            if file_path:
                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Export Successful", f"Plot has been saved to:\n{os.path.abspath(file_path)}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"An error occurred while exporting the plot:\n{str(e)}")

    def monomer_conversion(self, time, KP, KD, KT, INITIATOR_EFFICIENCY, INITIATOR_CONC0):
        pre = (2 * KP / KD) * np.sqrt((INITIATOR_EFFICIENCY * KD * INITIATOR_CONC0) / KT)
        exp_term = np.exp(-KD * time / 2) - 1
        return np.exp(pre * exp_term)

    def calc_molecular_weight(self, time, monomer_ratio, MONOMER_CONC0, MONOMER_MW, 
                             INITIATOR_CONC0, KD, INITIATOR_EFFICIENCY, 
                             CHAIN_TRANSFER_CONSTANT):
        monomer_conc = monomer_ratio * MONOMER_CONC0
        initiator_conc = INITIATOR_CONC0 * np.exp(-KD * time)
        
        consumed_initiator = INITIATOR_CONC0 - initiator_conc
        
        with np.errstate(divide='ignore', invalid='ignore'):
            mw = ((MONOMER_CONC0 - monomer_conc) /
                  (INITIATOR_EFFICIENCY * CHAIN_TRANSFER_CONSTANT * consumed_initiator)) * MONOMER_MW
        
        mw = np.nan_to_num(mw, nan=0.0)
        return mw
    
    def update_plot(self):
        try:
            # 获取当前参数
            temp = int(self.temp_var.get())
            monomer_conc = float(self.monomer_conc_var.get())
            initiator_conc = float(self.initiator_conc_var.get())
            termination_mode = self.termination_var.get()
            
            # 基于温度更新速率常数
            if temp == 50:
                kp = self.arrhenius_data[50]['KP']
                kt = self.arrhenius_data[50]['KT'] * termination_mode
                kd = self.arrhenius_data[50]['KD']
            else:
                R = 1.987
                T_K = temp + 273.15
                arr_data = self.arrhenius_data[70]
                kp = arr_data['KP_A'] * np.exp(-arr_data['KP_Ea'] / (R * T_K))
                kt = arr_data['KT_A'] * np.exp(-arr_data['KT_Ea'] / (R * T_K)) * termination_mode
                kd = arr_data['KD_A'] * np.exp(-arr_data['KD_Ea'] / (R * T_K))
            
            # 计算结果
            time_hours = np.arange(0, self.params['TIME_HOURS_MAX'], self.params['TIME_STEP'])
            time_seconds = time_hours * 3600

            monomer_ratio = self.monomer_conversion(time_seconds, kp, kd, kt, 
                                                   self.params['INITIATOR_EFFICIENCY'], 
                                                   initiator_conc)
            conversion = (1 - monomer_ratio) * 100
            molecular_weight = self.calc_molecular_weight(time_seconds, monomer_ratio, 
                                                         monomer_conc, self.params['MONOMER_MW'],
                                                         initiator_conc, kd,
                                                         self.params['INITIATOR_EFFICIENCY'], 
                                                         termination_mode)
            
            # 更新图表
            self.ax1.clear()
            self.ax2.clear()
            
            # 转化率图表
            self.ax1.plot(time_hours, conversion, 'r-', linewidth=2, label='Conversion')
            self.ax1.set_xlabel('Time (hours)')
            self.ax1.set_ylabel('Conversion (%)', color='red')
            self.ax1.tick_params(axis='y', labelcolor='red')
            self.ax1.set_ylim(0, 100)
            self.ax1.grid(True, alpha=0.3)
            self.ax1.legend(loc='lower right')
            
            # 分子量图表
            self.ax2.plot(time_hours, molecular_weight, 'b-', linewidth=2, label='Molecular Weight')
            self.ax2.set_xlabel('Time (hours)')
            self.ax2.set_ylabel('Molecular Weight (g/mol)', color='blue')
            self.ax2.tick_params(axis='y', labelcolor='blue')
            if max(molecular_weight) > 0:
                self.ax2.set_ylim(0, max(molecular_weight) * 1.1)
            else:
                self.ax2.set_ylim(0, 1000)
            self.ax2.grid(True, alpha=0.3)
            self.ax2.legend(loc='upper left')
            
            # 设置标题
            self.fig.suptitle(f'Polymerization Kinetics: {self.params["MONOMER"]} with {self.params["INITIATOR"]} at {temp} °C\n'
                             f'[M]₀ = {monomer_conc} M, [I]₀ = {initiator_conc} M, Termination mode = {termination_mode:.1f}')
            
            self.canvas.draw()
            self.update_parameters_display()
            
        except ValueError as e:
            messagebox.showerror("Input Error", "Please enter valid numerical values for all parameters.")

def main():
    # 设置matplotlib后端
    matplotlib.use('TkAgg')
    
    root = tk.Tk()
    app = PolymerizationKineticsApp(root)
    
    # 确保程序正确退出
    try:
        root.mainloop()
    except KeyboardInterrupt:
        app.on_closing()

if __name__ == "__main__":
    main()

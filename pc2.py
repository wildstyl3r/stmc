import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import sys

def find_column_by_prefix(df, prefix, file_label):
    matches = [col for col in df.columns if str(col).strip().startswith(prefix)]
    if len(matches) == 0:
        print(f"Ошибка: В файле [{file_label}] нет колонок, начинающихся с префикса '{prefix}'.")
        sys.exit(1)
    elif len(matches) > 1:
        print(f"Ошибка: В файле [{file_label}] найдено несколько колонок с префиксом '{prefix}': {matches}.")
        sys.exit(1)
    return matches[0]

def process_paschen_network(file_paths, ep_prefix, alpha_prefix, gamma_list):
    plt.figure(figsize=(10, 8))
    
    # Стили линий для разных значений gamma, чтобы отличать их внутри одного файла
    line_styles = ['-', '--', '-.', ':', (0, (3, 5, 1, 5))]
    
    # Цветовая палитра для разделения разных файлов
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    # Общая структура для сбора сводных данных в один CSV
    consolidated_export = {}

    print("=== Начало расчета кривых Пашена ===")
    
    for f_idx, path in enumerate(file_paths):
        if not os.path.exists(path):
            print(f"Ошибка: Файл не найден по пути '{path}'")
            sys.exit(1)
            
        # Извлекаем красивое имя файла для легенды и логов
        file_label = os.path.basename(path)
        print(f"\nОбработка файла: {file_label}")
        
        # Чтение данных
        df = pd.read_csv(path)
        
        # Динамический поиск колонок для текущего файла
        ep_col = find_column_by_prefix(df, ep_prefix, file_label)
        alpha_col = find_column_by_prefix(df, alpha_prefix, file_label)
        
        print(f"  -> Найдена колонка E/p: '{ep_col}'")
        print(f"  -> Найдена колонка alpha/p: '{alpha_col}'")
        
        # Конвертация в массивы
        ep = pd.to_numeric(df[ep_col], errors='coerce').values
        alpha_p = pd.to_numeric(df[alpha_col], errors='coerce').values
        
        # Очистка данных от некорректных значений и нулей
        valid_idx = (ep > 0) & (alpha_p > 0) & (~np.isnan(ep)) & (~np.isnan(alpha_p))
        ep = ep[valid_idx]
        alpha_p = alpha_p[valid_idx]
        
        if len(ep) == 0:
            print(f"  Ошибка: В файле '{file_label}' нет валидных числовых данных.")
            sys.exit(1)
            
        # Сортировка по возрастанию E/p для корректного построения линий графика
        sort_idx = np.argsort(ep)
        ep = ep[sort_idx]
        alpha_p = alpha_p[sort_idx]
        
        # Базовые колонки для экспорта
        consolidated_export[f'{file_label}_E_p'] = ep
        consolidated_export[f'{file_label}_alpha_p'] = alpha_p
        
        file_color = color_cycle[f_idx % len(color_cycle)]
        
        # Расчет для каждого значения gamma
        for g_idx, g in enumerate(sorted(gamma_list)):
            breakdown_criterion = np.log(1.0 + 1.0 / g)
            
            # Прямой расчет параметров пробоя
            pd_vals = breakdown_criterion / alpha_p
            v_break = ep * pd_vals
            
            # Сохраняем расчетные векторы в структуру экспорта
            consolidated_export[f'{file_label}_pd_(g={g})'] = pd_vals
            consolidated_export[f'{file_label}_Vb_(g={g})'] = v_break
            
            # Поиск минимума Пашена
            min_idx = np.argmin(v_break)
            print(f"    g = {g:<6} -> Минимум: pd = {pd_vals[min_idx]:.4f} Torr*cm, Vb = {v_break[min_idx]:.2f} V")
            
            # Настройка стиля линии
            l_style = line_styles[g_idx % len(line_styles)]
            
            # Отрисовка кривой
            plt.plot(
                pd_vals, 
                v_break, 
                marker='o', 
                linestyle=l_style, 
                color=file_color, 
                linewidth=1.8,
                label=f"{path.split('/')[-2].split("elasticCS-")[-1]} ($\gamma$={g})"
            )
            
    # Сохранение результатов
    # Так как массивы из разных файлов могут быть разной длины, 
    # оборачиваем в Series, чтобы pandas автоматически заполнил пустоты NaN
    export_series = {k: pd.Series(v) for k, v in consolidated_export.items()}
    output_df = pd.DataFrame(export_series)
    output_filename = "paschen_network_output.csv"
    output_df.to_csv(output_filename, index=False)
    print(f"\nГотово! Все расчетные данные экспортированы в '{output_filename}'")
    
    # Стилизация общего графика
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$p \cdot d$ (Torr$\cdot$cm)', fontsize=12)
    plt.ylabel('$V_b$ (V)', fontsize=12)
    # plt.title('Сводные кривые Пашена (Мульти-файл & Мульти-$\gamma$)', fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=9, loc='best', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.tight_layout()
    
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Использование: python paschen_multi_file_gamma.py <файлы_через_запятую> <E/p_префикс> <alpha/p_префикс> <gamma_через_запятую>")
        print("Пример: python paschen_multi_file_gamma.py he_data.csv,ar_data.csv E/p alpha 0.1,0.01")
        sys.exit(1)
        
    raw_files = sys.argv[1]
    ep_pref = sys.argv[2]
    alpha_pref = sys.argv[3]
    raw_gammas = sys.argv[4]
    
    # Парсинг списка файлов
    file_list = [f.strip() for f in raw_files.split(',') if f.strip()]
    if not file_list:
        print("Ошибка: Не указан ни один файл для обработки.")
        sys.exit(1)
        
    # Парсинг списка коэффициентов гамма
    gamma_list = []
    try:
        for item in raw_gammas.split(','):
            if not item.strip():
                continue
            g_val = float(item.strip())
            if g_val <= 0:
                print(f"Ошибка: Значение коэффициента гамма '{item}' должно быть строго больше нуля.")
                sys.exit(1)
            gamma_list.append(g_val)
    except ValueError:
        print(f"Ошибка: Не удалось распознать '{raw_gammas}' как список числовых значений float.")
        sys.exit(1)
        
    if not gamma_list:
        print("Ошибка: Не указано ни одно значение гамма.")
        sys.exit(1)
        
    process_paschen_network(file_list, ep_pref, alpha_pref, gamma_list)

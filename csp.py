import argparse
import sys
import matplotlib.pyplot as plt

# Валидные типы процессов LxCat
VALID_PROCESSES = {
    "ELASTIC", "IONIZATION", "EXCITATION", "EFFECTIVE", 
    "ROTATION", "DEEXCITATION", "ATTACHMENT",
    "ION_BACKSCATTER", "ION_ISOTROPIC"
}

def parse_and_filter_file(filepath, filter_str):
    """Парсит файл LxCat, извлекая отфильтрованные процессы, их описание и данные."""
    with open(filepath, 'r', encoding='utf-8') as f:
        # Сохраняем оригинальные строки (без \n), но не делаем strip заранее,
        # чтобы точно проверять пустые строки и структуру.
        lines = [line.rstrip('\r\n') for line in f]

    idx = 0
    num_lines = len(lines)

    while idx < num_lines:
        line_stripped = lines[idx].strip()
        
        # Проверяем точное совпадение с типом процесса
        if line_stripped in VALID_PROCESSES:
            process_type = line_stripped
            
            # Берем следующую строку как дополнительное описание (метаданные для легенды)
            process_details = ""
            if idx + 1 < num_lines:
                process_details = lines[idx + 1].strip()
                
            idx += 2
            
            # Ищем начало блока данных (разделитель из дефисов)
            while idx < num_lines:
                current_line = lines[idx].strip()
                if current_line.startswith('-----') and current_line.endswith('-----'):
                    idx += 1  # Переходим к самим данным
                    break
                idx += 1
            
            # Читаем координаты точек до первой пустой строки
            energies = []
            cross_sections = []
            
            while idx < num_lines and lines[idx].strip():
                parts = lines[idx].split()
                if len(parts) >= 2:
                    try:
                        energies.append(float(parts[0]))
                        cross_sections.append(float(parts[1]))
                    except ValueError:
                        # Игнорируем некорректные строки внутри блока данных
                        pass
                idx += 1
            
            # Проверяем фильтр (без учета регистра)
            if filter_str.lower() in process_type.lower():
                # Формируем полное имя для легенды: "Тип | Описание"
                full_description = f"{process_type} {process_details}".strip()
                yield full_description, energies, cross_sections
        else:
            idx += 1

def main():
    parser = argparse.ArgumentParser(
        description="Фильтрация и построение графиков сечений рассеяния электронов в формате LxCat/Bolsig."
    )
    parser.add_argument(
        "filter", 
        type=str, 
        help="Строка-фильтр для типа процесса (например, 'ION', 'EXC')."
    )
    parser.add_argument(
        "files", 
        type=str, 
        nargs="+", 
        help="Один или несколько текстовых файлов с данными."
    )
    
    args = parser.parse_args()
    
    plt.figure(figsize=(10, 6))
    data_plotted = False

    # Обрабатываем каждый переданный файл
    for filepath in args.files:
        try:
            for description, x, y in parse_and_filter_file(filepath, args.filter):
                if x and y:
                    # Добавляем имя файла и описание процесса в легенду
                    label = f"{filepath} ({description})"
                    plt.plot(x, y, label=label, marker='o', markersize=3)
                    data_plotted = True
        except FileNotFoundError:
            print(f"Ошибка: Файл не найден: {filepath}", file=sys.stderr)
        except Exception as e:
            print(f"Ошибка при обработке {filepath}: {e}", file=sys.stderr)

    # Отображаем график, если найдены соответствия
    if data_plotted:
        plt.xlabel("Energy (eV)")
        plt.ylabel("Cross Section (m²)")
        plt.yscale("log")
        plt.xscale("log")
        plt.title(f"Cross Sections Matching Filter: '{args.filter}'")
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.legend(loc="best")
        plt.tight_layout()
        plt.show()
    else:
        print(f"Не найдено сечений, соответствующих фильтру '{args.filter}'.")

if __name__ == "__main__":
    main()

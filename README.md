##hdreform
    **filename**
    **-o**, **--output** default=out - _optimal_
    **--tcf**, nargs=*, default= - _optimal_
    **-g**, **--gname**, nargs=*, default= - _optimal_

    Преобразовывает исходный hdf файл, полученный с помощью impulse. Высчитывает среднее по всем траекториям для указанных групп и типов корреляционных функций. 
    Группа tcf в атрибутах хранит group_size и names из исходного файла
    Полученный фал имеет следующую структуру: <group_name>/<tcf>/errors
                                                                /mean
                                              ...
                                              time

##hdfinfo
    **filename**, type=str
    **-a**, **--all** - _optimal_, статистика будет выдаваться по каждой группе дополнительно

    Проверяет файл на соответствие структуре файла получаемого с помощью impulse. Выдаёт статистику по файлу: число групп и траекторий, длина трейса и т.д.

##hdfextract
    filename
    --tcf 
    -o, --output, default=out, _optimal_
    -g, --gname,  nargs=*, default='' - _optimal_
    -m, --mean - _optimal_
    Извлекает определенную группу (или все) в файл .npy

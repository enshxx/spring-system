# spring-system

# Задача
Поставновку задачи смотрите в файле docs/problem.pdf

# Запуск 
Для запуска использовать из корневой папки ./build.sh; ./gauss_newton
Все параметры находятся в файле main.cpp
Целевые параметры по которым запускается симуляция - targetParams 
Начальные значения параметров - optParams, они будут подбираться в процессе 
работы алгоритма и в идеале должны сравняться с targetParams
Также есть параметры генерации экспериментальных точек для работы алгоритма
odeIntegrationStep, maxSteps, noiseLevel, 
paramStep - шаг для вычисления якобиана
Если модель расходится можно попробовать изменить начальное значения optParams

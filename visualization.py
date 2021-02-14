import matplotlib.pyplot as plt
from constants import delta_t, total_time, path


def plot_data(matr, layer, field, title='layer_field.png'):
    data = [getattr(matr[dt][layer], field) for dt in range(len(matr))]
    dt_list = [dt for dt in range(0, total_time, delta_t)]
    plt.figure(figsize=(16,9))
    plt.plot(dt_list, data, 'g')
    plt.xlabel('time')
    plt.ylabel(field)
    plt.show()
    try:
        plt.savefig(path + title)
    except:
        print('Please change path variable in constants to existing path')

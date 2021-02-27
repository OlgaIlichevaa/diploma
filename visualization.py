import matplotlib.pyplot as plt
from constants import delta_t, total_time, path


def plot_data(matr, layer, field, title='layer_field.png'):
    data = [getattr(matr[dt][layer], field) for dt in range(len(matr))]
    dt_list = [dt for dt in range(0, total_time, delta_t)]
    plt.figure(figsize=(16, 9))
    plt.plot(dt_list, data, 'g')
    plt.xlabel('time')
    plt.ylabel(field)
    #plt.show()
    try:
        plt.savefig(path + title)
    except:
        print('Please change path variable in constants to existing path')


def plot_all_layers(matr, n_layers, field, dt_list, title='all_layers.png'):
    layers_list = [layer for layer in range(n_layers)]
    plt.figure(figsize=(16, 9))
    for dt in dt_list:
        data_list = [getattr(matr[dt][layer], field) for layer in range(n_layers)]
        plt.plot(layers_list, data_list, label='dt = %s s' % (dt * delta_t))
    plt.xlabel('layer')
    plt.ylabel(field)
    plt.legend(loc='lower right')
    #plt.show()
    try:
        plt.savefig(path + title)
    except:
        print('Please change path variable in constants to existing path')

def plot_all_layers_no_zero_layer(matr, n_layers, field, dt_list, title='all_layers.png'):
    layers_list = [layer for layer in range(1, n_layers)]
    plt.figure(figsize=(16, 9))
    for dt in dt_list:
        data_list = [getattr(matr[dt][layer], field) for layer in range(1, n_layers)]
        plt.plot(layers_list, data_list, label='dt = %s s' % (dt * delta_t))
    plt.xlabel('layer')
    plt.ylabel(field)
    plt.legend(loc='lower right')
    #plt.show()
    try:
        plt.savefig(path + title)
    except:
        print('Please change path variable in constants to existing path')


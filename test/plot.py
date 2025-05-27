import os
import pandas as pd
import matplotlib.pyplot as plt

def grid_size_test(filename: str):
    df = pd.read_csv('data/' + filename)
    df['h'] = 1 / (df['n'] - 1)

    # First plot: x = n
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('n')
    ax1.set_ylabel('Time (s)')
    ax1.plot(df['n'], df['serial'], 'o-', color='tab:blue', label='Serial Time')
    ax1.plot(df['n'], df['omp'], 'o-', color='tab:orange', label='OMP Time')
    ax1.plot(df['n'], df['mpi'], 'o-', color='tab:green', label='MPI Time')
    ax1.plot(df['n'], df['hybrid'], 'o-', color='tab:purple', label='Hybrid Time')
    ax1.plot(df['n'], df['direct'], 'o-', color='tab:pink', label='Direct Time')
    ax1.set_xscale('log', base=2)
    ax1.set_yscale('log', base=2)
    ax1.legend(loc='upper center')
    ax2 = ax1.twinx()
    ax2.set_ylabel('L2 Error')
    ax2.plot(df['n'], df['l2_error'], 's-', color='tab:red', label='L2 Error')
    ax2.set_yscale('log', base=2)
    ax2.legend(loc='lower center')
    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
    plt.title('Timing and L2 Error vs Grid Size (n)')
    plt.tight_layout()
    plt.show()

    # Second plot: x = h
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('h = 1/(n-1)')
    ax1.set_ylabel('Time (s)')
    ax1.plot(df['h'], df['serial'], 'o-', color='tab:blue', label='Serial Time')
    ax1.plot(df['h'], df['omp'], 'o-', color='tab:orange', label='OMP Time')
    ax1.plot(df['h'], df['mpi'], 'o-', color='tab:green', label='MPI Time')
    ax1.plot(df['h'], df['hybrid'], 'o-', color='tab:purple', label='Hybrid Time')
    ax1.plot(df['h'], df['direct'], 'o-', color='tab:pink', label='Direct Time')
    ax1.set_xscale('log', base=2)
    ax1.set_yscale('log', base=2)
    ax1.legend(loc='upper center')
    ax2 = ax1.twinx()
    ax2.set_ylabel('L2 Error')
    ax2.plot(df['h'], df['l2_error'], 's-', color='tab:red', label='L2 Error')
    ax2.set_yscale('log', base=2)
    ax2.legend(loc='lower center')
    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
    plt.title('Timing and L2 Error vs Grid Spacing (h)')
    plt.tight_layout()
    plt.show()

grid_size_test("results_1.csv")
grid_size_test("results_2.csv")
grid_size_test("results_4.csv")
grid_size_test("results_8.csv")

def scalability_test():
    x = [1, 2, 4, 8]
    for i in [56, 64]:
        timings = []
        for filename in os.listdir('data'):
            if filename.endswith('.csv'):
                df = pd.read_csv('data/' + filename)
                df = df[df['n'] == i]
                timings.append(df['hybrid'].iloc[0])
        print(timings)
        plt.plot(x, timings, 'o-', label=f'n={i}')
        timings = []
    plt.xscale('log', base=2)
    plt.legend()
    plt.xlabel('Number of Processes')
    plt.ylabel('Time (s)')
    plt.title('Scalability Test')
    plt.tight_layout()
    plt.show()

scalability_test()
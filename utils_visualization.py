from matplotlib import pyplot as plt

def plot_images(imlist, cols=3, titlelist=None):
    n = imlist.shape[0]
    rows = int((n + cols - 1) / cols)
    for i, im in enumerate(imlist):
        plt.subplot(rows, cols, i + 1)
        plt.imshow(im)
        if titlelist is not None:
            plt.title(titlelist[i])
    plt.show()

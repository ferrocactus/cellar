import torch.nn as nn


class Autoencoder(nn.Module):
    """
    Module for producing a low dimensional representation for a set of vectors.

    Args:
        input_dim : Integer denoting the input dimensionality
        factor : Factor to reduce the width of successive layers in the encoder
                 Between 0 and 1 only (strictly)
        output_dim : Integer denoting the dimensionality of encoding
        non_linearity : Non -linearity to use between the layers. Default: None
                        If ``None``, then no non-linearity is added.
    """

    def __init__(self, input_dim, factor, output_dim, non_linearity=None):
        super(Autoencoder, self).__init__()
        assert 0 < factor < 1, "Factor value out of range"
        assert input_dim > output_dim, "Encoding dim cannot be greater than original dim"

        layers = []
        while int(factor * input_dim) > output_dim:
            layers.append((input_dim, int(factor * input_dim)))
            input_dim = int(input_dim * factor)
        layers.append((input_dim, output_dim))

        self.encoder = nn.Sequential()
        for i, l in enumerate(layers):
            self.encoder.add_module(
                'Linear_{}'.format(i), nn.Linear(l[0], l[1]))
            if i < len(layers) - 1 and non_linearity is not None:
                self.encoder.add_module('{}_{}'.format(non_linearity, i),
                                        getattr(nn, non_linearity)())

        self.decoder = nn.Sequential()
        for i, l in enumerate(reversed(layers)):
            self.decoder.add_module(
                'Linear_{}'.format(i), nn.Linear(l[1], l[0]))
            if i < len(layers) - 1 and non_linearity is not None:
                self.decoder.add_module('{}_{}'.format(non_linearity, i),
                                        getattr(nn, non_linearity)())

    def forward(self, input):
        """
        Function to perform the forward pass.

        Args:
            input : This is a tensor of inputs (batched)

        Returns:
            Reconstructed output
        """
        return self.decoder(self.encoder(input))

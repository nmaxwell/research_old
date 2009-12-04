
import pylab as p
X = [0,0.00390625,0.0078125,0.0117188,0.015625,0.0195312,0.0234375,0.0273438,0.03125,0.0351562,0.0390625,0.0429688,0.046875,0.0507812,0.0546875,0.0585938,0.0625,0.0664062,0.0703125,0.0742188,0.078125,0.0820312,0.0859375,0.0898438,0.09375,0.0976562,0.101562,0.105469,0.109375,0.113281,0.117188,0.121094,0.125,0.128906,0.132812,0.136719,0.140625,0.144531,0.148438,0.152344,0.15625,0.160156,0.164062,0.167969,0.171875,0.175781,0.179688,0.183594,0.1875,0.191406,0.195312,0.199219,0.203125,0.207031,0.210938,0.214844,0.21875,0.222656,0.226562,0.230469,0.234375,0.238281,0.242188,0.246094,0.25,0.253906,0.257812,0.261719,0.265625,0.269531,0.273438,0.277344,0.28125,0.285156,0.289062,0.292969,0.296875,0.300781,0.304688,0.308594,0.3125,0.316406,0.320312,0.324219,0.328125,0.332031,0.335938,0.339844,0.34375,0.347656,0.351562,0.355469,0.359375,0.363281,0.367188,0.371094,0.375,0.378906,0.382812,0.386719,0.390625,0.394531,0.398438,0.402344,0.40625,0.410156,0.414062,0.417969,0.421875,0.425781,0.429688,0.433594,0.4375,0.441406,0.445312,0.449219,0.453125,0.457031,0.460938,0.464844,0.46875,0.472656,0.476562,0.480469,0.484375,0.488281,0.492188,0.496094,0.5,0.503906,0.507812,0.511719,0.515625,0.519531,0.523438,0.527344,0.53125,0.535156,0.539062,0.542969,0.546875,0.550781,0.554688,0.558594,0.5625,0.566406,0.570312,0.574219,0.578125,0.582031,0.585938,0.589844,0.59375,0.597656,0.601562,0.605469,0.609375,0.613281,0.617188,0.621094,0.625,0.628906,0.632812,0.636719,0.640625,0.644531,0.648438,0.652344,0.65625,0.660156,0.664062,0.667969,0.671875,0.675781,0.679688,0.683594,0.6875,0.691406,0.695312,0.699219,0.703125,0.707031,0.710938,0.714844,0.71875,0.722656,0.726562,0.730469,0.734375,0.738281,0.742188,0.746094,0.75,0.753906,0.757812,0.761719,0.765625,0.769531,0.773438,0.777344,0.78125,0.785156,0.789062,0.792969,0.796875,0.800781,0.804688,0.808594,0.8125,0.816406,0.820312,0.824219,0.828125,0.832031,0.835938,0.839844,0.84375,0.847656,0.851562,0.855469,0.859375,0.863281,0.867188,0.871094,0.875,0.878906,0.882812,0.886719,0.890625,0.894531,0.898438,0.902344,0.90625,0.910156,0.914062,0.917969,0.921875,0.925781,0.929688,0.933594,0.9375,0.941406,0.945312,0.949219,0.953125,0.957031,0.960938,0.964844,0.96875,0.972656,0.976562,0.980469,0.984375,0.988281,0.992188,0.996094,1,1.00391,1.00781,1.01172,1.01562,1.01953,1.02344,1.02734,1.03125,1.03516,1.03906,1.04297,1.04688,1.05078,1.05469,1.05859,1.0625,1.06641,1.07031,1.07422,1.07812,1.08203,1.08594,1.08984,1.09375,1.09766,1.10156,1.10547,1.10938,1.11328,1.11719,1.12109,1.125,1.12891,1.13281,1.13672,1.14062,1.14453,1.14844,1.15234,1.15625,1.16016,1.16406,1.16797,1.17188,1.17578,1.17969,1.18359,1.1875,1.19141,1.19531,1.19922,1.20312,1.20703,1.21094,1.21484,1.21875,1.22266,1.22656,1.23047,1.23438,1.23828,1.24219,1.24609,1.25,1.25391,1.25781,1.26172,1.26562,1.26953,1.27344,1.27734,1.28125,1.28516,1.28906,1.29297,1.29688,1.30078,1.30469,1.30859,1.3125,1.31641,1.32031,1.32422,1.32812,1.33203,1.33594,1.33984,1.34375,1.34766,1.35156,1.35547,1.35938,1.36328,1.36719,1.37109,1.375,1.37891,1.38281,1.38672,1.39062,1.39453,1.39844,1.40234,1.40625,1.41016,1.41406,1.41797,1.42188,1.42578,1.42969,1.43359,1.4375,1.44141,1.44531,1.44922,1.45312,1.45703,1.46094,1.46484,1.46875,1.47266,1.47656,1.48047,1.48438,1.48828,1.49219,1.49609,1.5,1.50391,1.50781,1.51172,1.51562,1.51953,1.52344,1.52734,1.53125,1.53516,1.53906,1.54297,1.54688,1.55078,1.55469,1.55859,1.5625,1.56641,1.57031,1.57422,1.57812,1.58203,1.58594,1.58984,1.59375,1.59766,1.60156,1.60547,1.60938,1.61328,1.61719,1.62109,1.625,1.62891,1.63281,1.63672,1.64062,1.64453,1.64844,1.65234,1.65625,1.66016,1.66406,1.66797,1.67188,1.67578,1.67969,1.68359,1.6875,1.69141,1.69531,1.69922,1.70312,1.70703,1.71094,1.71484,1.71875,1.72266,1.72656,1.73047,1.73438,1.73828,1.74219,1.74609,1.75,1.75391,1.75781,1.76172,1.76562,1.76953,1.77344,1.77734,1.78125,1.78516,1.78906,1.79297,1.79688,1.80078,1.80469,1.80859,1.8125,1.81641,1.82031,1.82422,1.82812,1.83203,1.83594,1.83984,1.84375,1.84766,1.85156,1.85547,1.85938,1.86328,1.86719,1.87109,1.875,1.87891,1.88281,1.88672,1.89062,1.89453,1.89844,1.90234,1.90625,1.91016,1.91406,1.91797,1.92188,1.92578,1.92969,1.93359,1.9375,1.94141,1.94531,1.94922,1.95312,1.95703,1.96094,1.96484,1.96875,1.97266,1.97656,1.98047,1.98438,1.98828,1.99219,1.99609,2,2.00391,2.00781,2.01172,2.01562,2.01953,2.02344,2.02734,2.03125,2.03516,2.03906,2.04297,2.04688,2.05078,2.05469,2.05859,2.0625,2.06641,2.07031,2.07422,2.07812,2.08203,2.08594,2.08984,2.09375,2.09766,2.10156,2.10547,2.10938,2.11328,2.11719,2.12109,2.125,2.12891,2.13281,2.13672,2.14062,2.14453,2.14844,2.15234,2.15625,2.16016,2.16406,2.16797,2.17188,2.17578,2.17969,2.18359,2.1875,2.19141,2.19531,2.19922,2.20312,2.20703,2.21094,2.21484,2.21875,2.22266,2.22656,2.23047,2.23438,2.23828,2.24219,2.24609,2.25,2.25391,2.25781,2.26172,2.26562,2.26953,2.27344,2.27734,2.28125,2.28516,2.28906,2.29297,2.29688,2.30078,2.30469,2.30859,2.3125,2.31641,2.32031,2.32422,2.32812,2.33203,2.33594,2.33984,2.34375,2.34766,2.35156,2.35547,2.35938,2.36328,2.36719,2.37109,2.375,2.37891,2.38281,2.38672,2.39062,2.39453,2.39844,2.40234,2.40625,2.41016,2.41406,2.41797,2.42188,2.42578,2.42969,2.43359,2.4375,2.44141,2.44531,2.44922,2.45312,2.45703,2.46094,2.46484,2.46875,2.47266,2.47656,2.48047,2.48438,2.48828,2.49219,2.49609,2.5,2.50391,2.50781,2.51172,2.51562,2.51953,2.52344,2.52734,2.53125,2.53516,2.53906,2.54297,2.54688,2.55078,2.55469,2.55859,2.5625,2.56641,2.57031,2.57422,2.57812,2.58203,2.58594,2.58984,2.59375,2.59766,2.60156,2.60547,2.60938,2.61328,2.61719,2.62109,2.625,2.62891,2.63281,2.63672,2.64062,2.64453,2.64844,2.65234,2.65625,2.66016,2.66406,2.66797,2.67188,2.67578,2.67969,2.68359,2.6875,2.69141,2.69531,2.69922,2.70312,2.70703,2.71094,2.71484,2.71875,2.72266,2.72656,2.73047,2.73438,2.73828,2.74219,2.74609,2.75,2.75391,2.75781,2.76172,2.76562,2.76953,2.77344,2.77734,2.78125,2.78516,2.78906,2.79297,2.79688,2.80078,2.80469,2.80859,2.8125,2.81641,2.82031,2.82422,2.82812,2.83203,2.83594,2.83984,2.84375,2.84766,2.85156,2.85547,2.85938,2.86328,2.86719,2.87109,2.875,2.87891,2.88281,2.88672,2.89062,2.89453,2.89844,2.90234,2.90625,2.91016,2.91406,2.91797,2.92188,2.92578,2.92969,2.93359,2.9375,2.94141,2.94531,2.94922,2.95312,2.95703,2.96094,2.96484,2.96875,2.97266,2.97656,2.98047,2.98438,2.98828,2.99219,2.99609,3,3.00391,3.00781,3.01172,3.01562,3.01953,3.02344,3.02734,3.03125,3.03516,3.03906,3.04297,3.04688,3.05078,3.05469,3.05859,3.0625,3.06641,3.07031,3.07422,3.07812,3.08203,3.08594,3.08984,3.09375,3.09766,3.10156,3.10547,3.10938,3.11328,3.11719,3.12109,3.125,3.12891,3.13281,3.13672,3.14062,3.14453,3.14844,3.15234,3.15625,3.16016,3.16406,3.16797,3.17188,3.17578,3.17969,3.18359,3.1875,3.19141,3.19531,3.19922,3.20312,3.20703,3.21094,3.21484,3.21875,3.22266,3.22656,3.23047,3.23438,3.23828,3.24219,3.24609,3.25,3.25391,3.25781,3.26172,3.26562,3.26953,3.27344,3.27734,3.28125,3.28516,3.28906,3.29297,3.29688,3.30078,3.30469,3.30859,3.3125,3.31641,3.32031,3.32422,3.32812,3.33203,3.33594,3.33984,3.34375,3.34766,3.35156,3.35547,3.35938,3.36328,3.36719,3.37109,3.375,3.37891,3.38281,3.38672,3.39062,3.39453,3.39844,3.40234,3.40625,3.41016,3.41406,3.41797,3.42188,3.42578,3.42969,3.43359,3.4375,3.44141,3.44531,3.44922,3.45312,3.45703,3.46094,3.46484,3.46875,3.47266,3.47656,3.48047,3.48438,3.48828,3.49219,3.49609,3.5,3.50391,3.50781,3.51172,3.51562,3.51953,3.52344,3.52734,3.53125,3.53516,3.53906,3.54297,3.54688,3.55078,3.55469,3.55859,3.5625,3.56641,3.57031,3.57422,3.57812,3.58203,3.58594,3.58984,3.59375,3.59766,3.60156,3.60547,3.60938,3.61328,3.61719,3.62109,3.625,3.62891,3.63281,3.63672,3.64062,3.64453,3.64844,3.65234,3.65625,3.66016,3.66406,3.66797,3.67188,3.67578,3.67969,3.68359,3.6875,3.69141,3.69531,3.69922,3.70312,3.70703,3.71094,3.71484,3.71875,3.72266,3.72656,3.73047,3.73438,3.73828,3.74219,3.74609,3.75,3.75391,3.75781,3.76172,3.76562,3.76953,3.77344,3.77734,3.78125,3.78516,3.78906,3.79297,3.79688,3.80078,3.80469,3.80859,3.8125,3.81641,3.82031,3.82422,3.82812,3.83203,3.83594,3.83984,3.84375,3.84766,3.85156,3.85547,3.85938,3.86328,3.86719,3.87109,3.875,3.87891,3.88281,3.88672,3.89062,3.89453,3.89844,3.90234,3.90625,3.91016,3.91406,3.91797,3.92188,3.92578,3.92969,3.93359,3.9375,3.94141,3.94531,3.94922,3.95312,3.95703,3.96094,3.96484,3.96875,3.97266,3.97656,3.98047,3.98438,3.98828,3.99219,3.99609,];
Y = [3.11695e-14,6.21725e-15,9.32587e-15,4.44089e-15,9.76996e-15,6.66134e-15,8.88178e-16,2.22045e-14,1.77636e-14,3.4639e-14,2.57572e-14,9.76996e-15,1.95399e-14,8.88178e-15,1.42109e-14,0,1.77636e-14,3.55271e-14,5.15143e-14,3.55271e-14,4.44089e-14,1.59872e-14,2.13163e-14,3.90799e-14,1.77636e-15,1.59872e-14,1.77636e-14,2.30926e-14,3.73035e-14,2.30926e-14,2.30926e-14,2.4869e-14,5.32907e-15,2.66454e-14,2.30926e-14,2.66454e-14,4.08562e-14,2.84217e-14,3.01981e-14,4.44089e-14,4.26326e-14,1.95399e-14,5.86198e-14,8.70415e-14,4.08562e-14,2.30926e-14,1.42109e-14,5.15143e-14,9.41469e-14,1.95399e-14,3.01981e-14,8.88178e-15,7.99361e-15,6.4837e-14,1.59872e-14,3.10862e-14,4.26326e-14,1.42109e-14,3.01981e-14,2.22045e-15,8.43769e-15,9.54792e-15,1.26565e-14,4.11893e-14,1.38901e-14,5.66214e-15,7.99361e-15,1.9984e-15,8.88178e-15,2.08722e-14,1.77636e-15,2.22045e-14,3.81917e-14,2.66454e-15,2.30926e-14,5.32907e-14,2.22045e-14,8.88178e-15,1.77636e-15,3.55271e-15,3.01981e-14,7.28306e-14,1.01252e-13,9.9476e-14,6.21725e-14,1.59872e-14,2.30926e-14,7.99361e-14,9.23706e-14,9.59233e-14,9.05942e-14,1.42109e-14,9.41469e-14,1.03029e-13,1.01252e-13,9.23706e-14,1.77636e-14,9.05942e-14,1.04805e-13,1.03029e-13,9.05942e-14,2.84217e-14,9.23706e-14,1.19016e-13,9.41469e-14,8.34888e-14,3.55271e-15,3.19744e-14,8.52651e-14,8.52651e-14,7.81597e-14,7.81597e-14,8.88178e-15,5.32907e-14,7.99361e-14,5.15143e-14,4.44089e-14,1.59872e-14,5.50671e-14,3.10862e-14,5.68434e-14,2.93099e-14,1.06581e-14,1.15463e-14,2.62013e-14,1.28786e-14,5.55112e-15,1.12133e-14,2.10633e-15,1.14353e-14,1.39888e-14,5.77316e-15,1.82077e-14,1.90958e-14,2.4869e-14,1.33227e-14,4.70735e-14,6.21725e-15,2.22045e-14,1.59872e-14,8.88178e-14,4.79616e-14,5.9508e-14,6.03961e-14,5.32907e-15,8.70415e-14,8.88178e-14,6.57252e-14,6.03961e-14,8.88178e-15,2.66454e-14,1.01252e-13,9.59233e-14,9.76996e-14,9.59233e-14,1.06581e-14,8.88178e-14,1.06581e-13,1.01252e-13,8.88178e-14,3.01981e-14,4.08562e-14,4.9738e-14,1.95399e-14,1.29674e-13,2.43361e-13,2.38032e-13,1.27898e-13,5.32907e-15,3.90799e-14,3.55271e-14,6.75016e-14,3.37508e-14,3.55271e-14,5.32907e-15,7.81597e-14,1.54543e-13,1.82965e-13,8.79297e-14,4.61853e-14,2.84217e-14,3.55271e-15,3.28626e-14,3.73035e-14,1.77636e-15,4.08562e-14,2.70894e-14,1.11466e-13,1.02141e-14,2.90878e-14,1.33227e-15,6.32827e-15,1.88273e-15,1.22125e-14,1.26565e-14,2.15383e-14,5.77316e-15,9.23706e-14,6.88338e-14,6.03961e-14,8.88178e-15,4.44089e-15,2.4869e-14,4.26326e-14,9.76996e-15,1.68754e-14,1.08358e-13,1.84741e-13,1.9007e-13,1.06581e-13,3.90799e-14,4.79616e-14,2.66454e-14,2.30926e-14,3.01981e-14,4.79616e-14,2.30926e-14,1.29674e-13,2.43361e-13,2.43361e-13,1.22569e-13,8.88178e-15,5.15143e-14,4.26326e-14,4.61853e-14,4.26326e-14,4.44089e-14,1.77636e-14,1.35003e-13,2.38032e-13,2.54019e-13,1.19016e-13,2.84217e-14,4.9738e-14,1.59872e-14,5.86198e-14,2.66454e-14,3.90799e-14,2.4869e-14,8.34888e-14,1.47438e-13,1.75859e-13,9.32587e-14,2.39808e-14,7.10543e-15,1.06581e-14,3.10862e-14,9.76996e-15,3.4639e-14,4.52971e-14,8.52651e-14,6.66134e-14,1.19904e-14,3.10862e-15,1.68754e-14,3.67484e-14,1.85207e-14,1.02141e-14,1.02141e-14,2.02061e-14,4.44089e-15,6.57252e-14,7.01661e-14,1.77636e-15,3.28626e-14,1.86517e-14,2.66454e-15,4.79616e-14,4.35207e-14,1.06581e-14,8.61533e-14,1.95399e-13,1.72307e-13,9.23706e-14,3.37508e-14,5.15143e-14,2.4869e-14,2.66454e-14,3.55271e-14,3.73035e-14,0,1.20792e-13,2.32703e-13,2.41585e-13,1.27898e-13,1.95399e-14,4.08562e-14,4.44089e-14,4.44089e-14,4.79616e-14,5.32907e-14,1.42109e-14,1.27898e-13,2.36255e-13,2.46914e-13,1.15463e-13,2.13163e-14,5.32907e-14,1.95399e-14,4.79616e-14,3.19744e-14,6.21725e-14,3.19744e-14,8.88178e-14,1.7053e-13,1.98952e-13,1.01252e-13,2.30926e-14,3.64153e-14,1.24345e-14,1.06581e-14,1.5099e-14,3.37508e-14,4.70735e-14,8.34888e-14,6.4837e-14,1.33227e-14,4.66294e-15,9.99201e-15,3.51941e-14,2.00595e-14,1.26565e-14,8.88178e-15,4.70735e-14,2.57572e-14,4.66294e-14,1.2168e-13,8.61533e-14,4.61853e-14,2.04281e-14,1.44773e-13,4.21885e-13,1.85629e-13,6.12843e-14,8.08242e-14,8.88178e-15,5.50671e-14,5.32907e-15,8.88178e-14,3.55271e-14,2.84217e-14,2.25597e-13,7.10543e-13,2.93099e-13,3.55271e-14,6.57252e-14,3.19744e-14,5.32907e-14,3.37508e-14,5.15143e-14,5.32907e-14,2.30926e-14,3.05533e-13,8.10019e-13,2.98428e-13,2.13163e-14,5.32907e-14,6.03961e-14,3.01981e-14,6.57252e-14,3.73035e-14,2.98428e-13,7.17648e-13,2.30926e-13,5.32907e-15,2.66454e-14,4.79616e-14,1.95399e-14,7.81597e-14,2.66454e-14,9.59233e-14,6.12843e-14,1.69642e-13,4.02345e-13,1.22569e-13,2.66454e-14,9.76996e-15,4.44089e-14,3.24185e-14,5.4623e-14,4.44089e-14,6.41709e-14,5.01821e-14,1.90958e-14,1.41716e-14,3.56382e-14,7.77156e-14,6.77236e-14,2.70894e-14,2.39808e-14,4.84057e-14,6.92779e-14,1.15463e-14,1.86517e-14,1.15463e-13,4.16556e-13,2.08722e-13,6.21725e-14,8.70415e-14,1.77636e-14,4.44089e-14,8.88178e-15,6.57252e-14,2.66454e-14,3.55271e-15,2.36255e-13,7.17648e-13,2.6823e-13,3.73035e-14,8.70415e-14,2.84217e-14,5.15143e-14,4.26326e-14,2.4869e-14,3.03757e-13,8.08242e-13,3.03757e-13,2.66454e-14,5.86198e-14,4.26326e-14,2.84217e-14,5.86198e-14,2.84217e-14,7.28306e-14,5.86198e-14,2.80664e-13,6.98108e-13,2.36255e-13,3.55271e-15,1.95399e-14,5.50671e-14,1.77636e-15,8.52651e-14,8.88178e-15,6.66134e-14,5.32907e-14,1.75859e-13,3.97904e-13,1.03029e-13,1.06581e-14,7.10543e-15,1.33227e-14,3.01981e-14,5.59552e-14,7.50511e-14,7.99361e-14,5.21805e-14,2.15383e-14,1.27437e-14,3.40838e-14,7.92699e-14,6.61693e-14,2.88658e-14,2.26485e-14,2.13163e-14,9.68114e-14,4.08562e-14,1.68754e-14,1.18128e-13,4.20997e-13,2.15827e-13,5.77316e-14,7.63833e-14,7.10543e-15,4.9738e-14,3.01981e-14,3.90799e-14,1.77636e-15,1.77636e-14,2.57572e-13,7.35412e-13,2.82441e-13,3.19744e-14,7.99361e-14,8.88178e-15,5.32907e-14,4.08562e-14,5.15143e-14,5.50671e-14,3.01981e-14,3.03757e-13,8.10019e-13,2.98428e-13,1.95399e-14,5.15143e-14,4.44089e-14,2.84217e-14,8.17124e-14,5.15143e-14,2.84217e-13,7.19425e-13,2.43361e-13,1.77636e-14,5.32907e-15,4.61853e-14,1.77636e-14,7.10543e-14,1.24345e-14,5.24025e-14,4.08562e-14,2.0961e-13,4.25437e-13,8.43769e-14,1.59872e-14,2.30926e-14,7.99361e-14,2.17604e-14,2.62013e-14,7.41629e-14,8.14904e-14,4.9738e-14,2.29816e-14,1.11979e-14,3.25295e-14,8.08242e-14,6.4615e-14,3.01981e-14,7.81597e-14,2.30926e-14,6.66134e-14,0,1.24345e-14,8.79297e-14,4.15668e-13,2.10498e-13,5.68434e-14,7.54952e-14,5.32907e-15,4.79616e-14,1.77636e-15,5.32907e-14,1.77636e-15,8.88178e-15,2.6823e-13,7.14095e-13,2.8777e-13,4.26326e-14,7.4607e-14,2.84217e-14,5.15143e-14,4.44089e-14,2.30926e-14,3.00204e-13,8.08242e-13,3.03757e-13,2.66454e-14,5.50671e-14,4.61853e-14,3.19744e-14,5.32907e-14,3.55271e-14,6.39488e-14,5.86198e-14,2.94875e-13,7.10543e-13,2.39808e-13,3.01981e-14,1.06581e-14,2.4869e-14,3.55271e-15,7.63833e-14,1.24345e-14,6.03961e-14,4.08562e-14,1.83853e-13,3.89022e-13,8.17124e-14,1.06581e-14,1.5099e-14,7.81597e-14,2.04281e-14,5.63993e-14,1.55431e-14,5.4623e-14,6.23945e-14,2.45359e-14,8.77082e-15,3.03091e-14,8.23785e-14,6.32827e-14,3.15303e-14,7.63833e-14,2.44249e-14,6.4837e-14,3.37508e-14,3.90799e-14,1.22569e-13,4.16556e-13,1.88294e-13,5.50671e-14,8.17124e-14,3.55271e-15,5.15143e-14,1.06581e-14,5.68434e-14,5.32907e-15,1.42109e-14,2.45137e-13,7.2653e-13,2.91323e-13,3.01981e-14,6.57252e-14,1.95399e-14,5.68434e-14,3.01981e-14,4.9738e-14,6.03961e-14,2.84217e-14,3.01981e-13,8.08242e-13,2.98428e-13,3.01981e-14,4.26326e-14,4.9738e-14,3.55271e-14,6.92779e-14,5.50671e-14,2.84217e-13,7.05214e-13,2.43361e-13,5.32907e-15,3.01981e-14,5.86198e-14,5.32907e-15,5.32907e-14,5.32907e-15,9.50351e-14,3.37508e-14,2.18492e-13,4.10338e-13,1.19016e-13,1.77636e-15,1.68754e-14,8.70415e-14,2.26485e-14,7.19425e-14,4.61853e-14,5.70655e-14,6.10623e-14,2.60902e-14,7.34296e-15,5.71765e-14,2.68674e-14,6.17284e-14,6.17284e-14,4.66294e-14,2.57572e-14,6.12843e-14,6.30607e-14,8.08242e-14,2.13163e-14,1.58096e-13,3.33067e-13,3.54383e-13,4.24549e-13,5.31131e-13,5.20473e-13,4.583e-13,5.68434e-13,6.19949e-13,5.52447e-13,6.80345e-13,1.09246e-12,7.49623e-13,5.63105e-13,6.89226e-13,7.31859e-13,6.55476e-13,6.43041e-13,7.35412e-13,7.51399e-13,6.35936e-13,8.10019e-13,1.23812e-12,8.0469e-13,6.30607e-13,7.49623e-13,7.19425e-13,6.18172e-13,6.51923e-13,6.89226e-13,6.50147e-13,5.50671e-13,7.05214e-13,1.04272e-12,6.53699e-13,5.09814e-13,5.80869e-13,4.72511e-13,3.8014e-13,4.66294e-13,4.73399e-13,3.90799e-13,2.60236e-13,4.30767e-13,5.25802e-13,2.86882e-13,2.56684e-13,2.77112e-13,1.01696e-13,1.25233e-13,1.37668e-13,1.2701e-13,2.93099e-14,2.37416e-14,1.88738e-15,8.4821e-14,1.54543e-13,7.63833e-14,1.58984e-13,2.74891e-13,2.78888e-13,3.10862e-13,5.62217e-13,3.94351e-13,3.14415e-13,3.84581e-13,5.11591e-13,4.84945e-13,4.15668e-13,5.20473e-13,5.98632e-13,5.18696e-13,6.53699e-13,1.04095e-12,6.98108e-13,5.52447e-13,6.43041e-13,7.12319e-13,6.35936e-13,6.25278e-13,7.17648e-13,7.49623e-13,6.35936e-13,8.02913e-13,1.24345e-12,8.08242e-13,6.34159e-13,7.47846e-13,7.40741e-13,6.39488e-13,6.66134e-13,7.38964e-13,6.98108e-13,5.63105e-13,7.47846e-13,1.09956e-12,6.8745e-13,5.77316e-13,6.4837e-13,5.66658e-13,4.65405e-13,5.31131e-13,5.63105e-13,4.42313e-13,3.54383e-13,4.69846e-13,6.26166e-13,3.52607e-13,2.98428e-13,3.52607e-13,2.12275e-13,1.6076e-13,2.31815e-13,1.86073e-13,9.17044e-14,1.90958e-14,7.63833e-14,2.48712e-14,6.51701e-14,1.57652e-14,4.4631e-14,1.73195e-13,1.90514e-13,1.2923e-13,1.9984e-13,3.00204e-13,3.09974e-13,3.23297e-13,6.2883e-13,5.01821e-13,3.39284e-13,4.1922e-13,5.29354e-13,4.88498e-13,4.61853e-13,5.57776e-13,6.44818e-13,5.73763e-13,6.85674e-13,1.08891e-12,7.54952e-13,5.63105e-13,6.89226e-13,7.33635e-13,6.57252e-13,6.37712e-13,7.42517e-13,7.4607e-13,6.35936e-13,8.06466e-13,1.23634e-12,7.97584e-13,6.32383e-13,7.42517e-13,7.19425e-13,6.25278e-13,6.37712e-13,7.14095e-13,6.59028e-13,5.40012e-13,7.22977e-13,1.00897e-12,6.27054e-13,4.83169e-13,5.86198e-13,4.84945e-13,3.96128e-13,4.52083e-13,4.85834e-13,3.7037e-13,2.71783e-13,4.20997e-13,5.31131e-13,2.50466e-13,2.22045e-13,2.55351e-13,1.5854e-13,9.37028e-14,1.63203e-13,1.5854e-13,2.9754e-14,9.00036e-14,1.53211e-14,1.13687e-13,1.5854e-13,7.90479e-14,1.41664e-13,2.85993e-13,2.39808e-13,3.13527e-13,5.60441e-13,3.93463e-13,3.09086e-13,3.57936e-13,5.13367e-13,4.49418e-13,4.29878e-13,4.93827e-13,6.07514e-13,5.32907e-13,6.21725e-13,1.02851e-12,7.14095e-13,5.29354e-13,6.44818e-13,7.17648e-13,6.41265e-13,6.18172e-13,7.19425e-13,7.42517e-13,6.2883e-13,8.0469e-13,1.23279e-12,8.13571e-13,6.37712e-13,7.51399e-13,7.35412e-13,6.44818e-13,6.44818e-13,7.28306e-13,6.8745e-13,5.61329e-13,7.51399e-13,1.08713e-12,6.91003e-13,5.3646e-13,6.27054e-13,5.50671e-13,4.65405e-13,5.04485e-13,5.70211e-13,4.23661e-13,3.68594e-13,4.64517e-13,6.35936e-13,3.4639e-13,2.88658e-13,3.31291e-13,2.22933e-13,1.21236e-13,1.74527e-13,1.6076e-13,6.63913e-14,1.95399e-14,7.71605e-14,2.78449e-14,6.80567e-14,1.88738e-14,4.15223e-14,1.95399e-13,1.93623e-13,1.60316e-13,2.11386e-13,3.25073e-13,3.12639e-13,3.5083e-13,6.5814e-13,4.72511e-13,3.6593e-13,4.44977e-13,5.48894e-13,4.93827e-13,4.9738e-13,5.24025e-13,6.39488e-13,5.56e-13,6.59028e-13,1.08891e-12,7.42517e-13,5.59552e-13,6.80345e-13,7.22977e-13,6.4837e-13,6.46594e-13,7.31859e-13,7.53175e-13,6.35936e-13,8.08242e-13,1.23634e-12,7.95808e-13,6.37712e-13,7.33635e-13,7.28306e-13,6.32383e-13,6.46594e-13,7.14095e-13,6.57252e-13,5.3646e-13,6.96332e-13,1.03739e-12,6.16396e-13,5.45342e-13,5.84421e-13,5.27578e-13,4.49418e-13,4.67182e-13,4.64517e-13,3.82805e-13,2.93099e-13,4.11227e-13,5.40901e-13,2.72671e-13,2.40696e-13,2.66009e-13,1.19904e-13,7.90479e-14,1.45217e-13,1.29674e-13,1.16573e-14,2.72396e-14,1.4877e-14,1.10578e-13,1.61648e-13,9.63674e-14,1.63425e-13,2.83773e-13,2.42473e-13,2.73559e-13,5.25802e-13,3.88134e-13,2.79776e-13,3.59712e-13,4.57412e-13,4.55636e-13,4.05009e-13,4.88498e-13,5.68434e-13,4.92051e-13,6.05738e-13,1.02496e-12,7.05214e-13,5.34683e-13,6.66134e-13,7.17648e-13,6.55476e-13,6.12843e-13,7.24754e-13,7.38964e-13,6.37712e-13,7.99361e-13,1.23812e-12,8.10019e-13,6.35936e-13,7.54952e-13,7.4607e-13,6.53699e-13,6.51923e-13,7.42517e-13,6.89226e-13,5.61329e-13,7.33635e-13,1.07114e-12,6.91003e-13,5.3646e-13,6.25278e-13,5.50671e-13,4.3876e-13,5.20473e-13,5.25802e-13,4.36984e-13,3.33067e-13,4.84057e-13,6.22613e-13,3.33067e-13,2.80664e-13,3.09086e-13,2.01617e-13,1.43441e-13,1.89182e-13,1.61204e-13,7.01661e-14,1.64313e-14,8.82627e-14,];
p.plot(X,Y) 
p.show()

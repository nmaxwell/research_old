
import pylab as p
X = [0,0.00390625,0.0078125,0.0117188,0.015625,0.0195312,0.0234375,0.0273438,0.03125,0.0351562,0.0390625,0.0429688,0.046875,0.0507812,0.0546875,0.0585938,0.0625,0.0664062,0.0703125,0.0742188,0.078125,0.0820312,0.0859375,0.0898438,0.09375,0.0976562,0.101562,0.105469,0.109375,0.113281,0.117188,0.121094,0.125,0.128906,0.132812,0.136719,0.140625,0.144531,0.148438,0.152344,0.15625,0.160156,0.164062,0.167969,0.171875,0.175781,0.179688,0.183594,0.1875,0.191406,0.195312,0.199219,0.203125,0.207031,0.210938,0.214844,0.21875,0.222656,0.226562,0.230469,0.234375,0.238281,0.242188,0.246094,0.25,0.253906,0.257812,0.261719,0.265625,0.269531,0.273438,0.277344,0.28125,0.285156,0.289062,0.292969,0.296875,0.300781,0.304688,0.308594,0.3125,0.316406,0.320312,0.324219,0.328125,0.332031,0.335938,0.339844,0.34375,0.347656,0.351562,0.355469,0.359375,0.363281,0.367188,0.371094,0.375,0.378906,0.382812,0.386719,0.390625,0.394531,0.398438,0.402344,0.40625,0.410156,0.414062,0.417969,0.421875,0.425781,0.429688,0.433594,0.4375,0.441406,0.445312,0.449219,0.453125,0.457031,0.460938,0.464844,0.46875,0.472656,0.476562,0.480469,0.484375,0.488281,0.492188,0.496094,0.5,0.503906,0.507812,0.511719,0.515625,0.519531,0.523438,0.527344,0.53125,0.535156,0.539062,0.542969,0.546875,0.550781,0.554688,0.558594,0.5625,0.566406,0.570312,0.574219,0.578125,0.582031,0.585938,0.589844,0.59375,0.597656,0.601562,0.605469,0.609375,0.613281,0.617188,0.621094,0.625,0.628906,0.632812,0.636719,0.640625,0.644531,0.648438,0.652344,0.65625,0.660156,0.664062,0.667969,0.671875,0.675781,0.679688,0.683594,0.6875,0.691406,0.695312,0.699219,0.703125,0.707031,0.710938,0.714844,0.71875,0.722656,0.726562,0.730469,0.734375,0.738281,0.742188,0.746094,0.75,0.753906,0.757812,0.761719,0.765625,0.769531,0.773438,0.777344,0.78125,0.785156,0.789062,0.792969,0.796875,0.800781,0.804688,0.808594,0.8125,0.816406,0.820312,0.824219,0.828125,0.832031,0.835938,0.839844,0.84375,0.847656,0.851562,0.855469,0.859375,0.863281,0.867188,0.871094,0.875,0.878906,0.882812,0.886719,0.890625,0.894531,0.898438,0.902344,0.90625,0.910156,0.914062,0.917969,0.921875,0.925781,0.929688,0.933594,0.9375,0.941406,0.945312,0.949219,0.953125,0.957031,0.960938,0.964844,0.96875,0.972656,0.976562,0.980469,0.984375,0.988281,0.992188,0.996094,1,1.00391,1.00781,1.01172,1.01562,1.01953,1.02344,1.02734,1.03125,1.03516,1.03906,1.04297,1.04688,1.05078,1.05469,1.05859,1.0625,1.06641,1.07031,1.07422,1.07812,1.08203,1.08594,1.08984,1.09375,1.09766,1.10156,1.10547,1.10938,1.11328,1.11719,1.12109,1.125,1.12891,1.13281,1.13672,1.14062,1.14453,1.14844,1.15234,1.15625,1.16016,1.16406,1.16797,1.17188,1.17578,1.17969,1.18359,1.1875,1.19141,1.19531,1.19922,1.20312,1.20703,1.21094,1.21484,1.21875,1.22266,1.22656,1.23047,1.23438,1.23828,1.24219,1.24609,1.25,1.25391,1.25781,1.26172,1.26562,1.26953,1.27344,1.27734,1.28125,1.28516,1.28906,1.29297,1.29688,1.30078,1.30469,1.30859,1.3125,1.31641,1.32031,1.32422,1.32812,1.33203,1.33594,1.33984,1.34375,1.34766,1.35156,1.35547,1.35938,1.36328,1.36719,1.37109,1.375,1.37891,1.38281,1.38672,1.39062,1.39453,1.39844,1.40234,1.40625,1.41016,1.41406,1.41797,1.42188,1.42578,1.42969,1.43359,1.4375,1.44141,1.44531,1.44922,1.45312,1.45703,1.46094,1.46484,1.46875,1.47266,1.47656,1.48047,1.48438,1.48828,1.49219,1.49609,1.5,1.50391,1.50781,1.51172,1.51562,1.51953,1.52344,1.52734,1.53125,1.53516,1.53906,1.54297,1.54688,1.55078,1.55469,1.55859,1.5625,1.56641,1.57031,1.57422,1.57812,1.58203,1.58594,1.58984,1.59375,1.59766,1.60156,1.60547,1.60938,1.61328,1.61719,1.62109,1.625,1.62891,1.63281,1.63672,1.64062,1.64453,1.64844,1.65234,1.65625,1.66016,1.66406,1.66797,1.67188,1.67578,1.67969,1.68359,1.6875,1.69141,1.69531,1.69922,1.70312,1.70703,1.71094,1.71484,1.71875,1.72266,1.72656,1.73047,1.73438,1.73828,1.74219,1.74609,1.75,1.75391,1.75781,1.76172,1.76562,1.76953,1.77344,1.77734,1.78125,1.78516,1.78906,1.79297,1.79688,1.80078,1.80469,1.80859,1.8125,1.81641,1.82031,1.82422,1.82812,1.83203,1.83594,1.83984,1.84375,1.84766,1.85156,1.85547,1.85938,1.86328,1.86719,1.87109,1.875,1.87891,1.88281,1.88672,1.89062,1.89453,1.89844,1.90234,1.90625,1.91016,1.91406,1.91797,1.92188,1.92578,1.92969,1.93359,1.9375,1.94141,1.94531,1.94922,1.95312,1.95703,1.96094,1.96484,1.96875,1.97266,1.97656,1.98047,1.98438,1.98828,1.99219,1.99609,2,2.00391,2.00781,2.01172,2.01562,2.01953,2.02344,2.02734,2.03125,2.03516,2.03906,2.04297,2.04688,2.05078,2.05469,2.05859,2.0625,2.06641,2.07031,2.07422,2.07812,2.08203,2.08594,2.08984,2.09375,2.09766,2.10156,2.10547,2.10938,2.11328,2.11719,2.12109,2.125,2.12891,2.13281,2.13672,2.14062,2.14453,2.14844,2.15234,2.15625,2.16016,2.16406,2.16797,2.17188,2.17578,2.17969,2.18359,2.1875,2.19141,2.19531,2.19922,2.20312,2.20703,2.21094,2.21484,2.21875,2.22266,2.22656,2.23047,2.23438,2.23828,2.24219,2.24609,2.25,2.25391,2.25781,2.26172,2.26562,2.26953,2.27344,2.27734,2.28125,2.28516,2.28906,2.29297,2.29688,2.30078,2.30469,2.30859,2.3125,2.31641,2.32031,2.32422,2.32812,2.33203,2.33594,2.33984,2.34375,2.34766,2.35156,2.35547,2.35938,2.36328,2.36719,2.37109,2.375,2.37891,2.38281,2.38672,2.39062,2.39453,2.39844,2.40234,2.40625,2.41016,2.41406,2.41797,2.42188,2.42578,2.42969,2.43359,2.4375,2.44141,2.44531,2.44922,2.45312,2.45703,2.46094,2.46484,2.46875,2.47266,2.47656,2.48047,2.48438,2.48828,2.49219,2.49609,2.5,2.50391,2.50781,2.51172,2.51562,2.51953,2.52344,2.52734,2.53125,2.53516,2.53906,2.54297,2.54688,2.55078,2.55469,2.55859,2.5625,2.56641,2.57031,2.57422,2.57812,2.58203,2.58594,2.58984,2.59375,2.59766,2.60156,2.60547,2.60938,2.61328,2.61719,2.62109,2.625,2.62891,2.63281,2.63672,2.64062,2.64453,2.64844,2.65234,2.65625,2.66016,2.66406,2.66797,2.67188,2.67578,2.67969,2.68359,2.6875,2.69141,2.69531,2.69922,2.70312,2.70703,2.71094,2.71484,2.71875,2.72266,2.72656,2.73047,2.73438,2.73828,2.74219,2.74609,2.75,2.75391,2.75781,2.76172,2.76562,2.76953,2.77344,2.77734,2.78125,2.78516,2.78906,2.79297,2.79688,2.80078,2.80469,2.80859,2.8125,2.81641,2.82031,2.82422,2.82812,2.83203,2.83594,2.83984,2.84375,2.84766,2.85156,2.85547,2.85938,2.86328,2.86719,2.87109,2.875,2.87891,2.88281,2.88672,2.89062,2.89453,2.89844,2.90234,2.90625,2.91016,2.91406,2.91797,2.92188,2.92578,2.92969,2.93359,2.9375,2.94141,2.94531,2.94922,2.95312,2.95703,2.96094,2.96484,2.96875,2.97266,2.97656,2.98047,2.98438,2.98828,2.99219,2.99609,3,3.00391,3.00781,3.01172,3.01562,3.01953,3.02344,3.02734,3.03125,3.03516,3.03906,3.04297,3.04688,3.05078,3.05469,3.05859,3.0625,3.06641,3.07031,3.07422,3.07812,3.08203,3.08594,3.08984,3.09375,3.09766,3.10156,3.10547,3.10938,3.11328,3.11719,3.12109,3.125,3.12891,3.13281,3.13672,3.14062,3.14453,3.14844,3.15234,3.15625,3.16016,3.16406,3.16797,3.17188,3.17578,3.17969,3.18359,3.1875,3.19141,3.19531,3.19922,3.20312,3.20703,3.21094,3.21484,3.21875,3.22266,3.22656,3.23047,3.23438,3.23828,3.24219,3.24609,3.25,3.25391,3.25781,3.26172,3.26562,3.26953,3.27344,3.27734,3.28125,3.28516,3.28906,3.29297,3.29688,3.30078,3.30469,3.30859,3.3125,3.31641,3.32031,3.32422,3.32812,3.33203,3.33594,3.33984,3.34375,3.34766,3.35156,3.35547,3.35938,3.36328,3.36719,3.37109,3.375,3.37891,3.38281,3.38672,3.39062,3.39453,3.39844,3.40234,3.40625,3.41016,3.41406,3.41797,3.42188,3.42578,3.42969,3.43359,3.4375,3.44141,3.44531,3.44922,3.45312,3.45703,3.46094,3.46484,3.46875,3.47266,3.47656,3.48047,3.48438,3.48828,3.49219,3.49609,3.5,3.50391,3.50781,3.51172,3.51562,3.51953,3.52344,3.52734,3.53125,3.53516,3.53906,3.54297,3.54688,3.55078,3.55469,3.55859,3.5625,3.56641,3.57031,3.57422,3.57812,3.58203,3.58594,3.58984,3.59375,3.59766,3.60156,3.60547,3.60938,3.61328,3.61719,3.62109,3.625,3.62891,3.63281,3.63672,3.64062,3.64453,3.64844,3.65234,3.65625,3.66016,3.66406,3.66797,3.67188,3.67578,3.67969,3.68359,3.6875,3.69141,3.69531,3.69922,3.70312,3.70703,3.71094,3.71484,3.71875,3.72266,3.72656,3.73047,3.73438,3.73828,3.74219,3.74609,3.75,3.75391,3.75781,3.76172,3.76562,3.76953,3.77344,3.77734,3.78125,3.78516,3.78906,3.79297,3.79688,3.80078,3.80469,3.80859,3.8125,3.81641,3.82031,3.82422,3.82812,3.83203,3.83594,3.83984,3.84375,3.84766,3.85156,3.85547,3.85938,3.86328,3.86719,3.87109,3.875,3.87891,3.88281,3.88672,3.89062,3.89453,3.89844,3.90234,3.90625,3.91016,3.91406,3.91797,3.92188,3.92578,3.92969,3.93359,3.9375,3.94141,3.94531,3.94922,3.95312,3.95703,3.96094,3.96484,3.96875,3.97266,3.97656,3.98047,3.98438,3.98828,3.99219,3.99609,];
Y = [0,1.21978e-17,4.21118e-17,4.0145e-17,-2.00402e-17,-2.069e-17,-4.38866e-17,2.33839e-17,-2.35032e-17,-7.53426e-19,1.0823e-17,3.73623e-17,-1.41645e-17,1.82422e-17,1.55711e-17,-4.56691e-18,2.66876e-17,-9.46362e-17,-3.70458e-17,9.5259e-17,3.21077e-17,-9.97763e-17,1.2289e-16,6.18259e-17,2.31334e-17,-4.46214e-17,-8.41922e-17,4.19017e-17,1.09534e-17,-5.86181e-17,1.11039e-16,5.24545e-17,-2.02227e-21,-5.68105e-17,-1.08014e-16,4.85013e-17,7.7129e-18,-4.9674e-17,1.33534e-16,4.03378e-17,-2.52538e-17,1.60534e-16,-1.19872e-16,-5.26514e-18,1.99134e-16,-8.26012e-17,5.34017e-17,-1.07126e-16,-2.26087e-18,1.44369e-16,-9.32618e-17,1.9828e-17,-1.64896e-16,1.07006e-17,4.24698e-17,-5.1607e-17,-2.3362e-17,4.63813e-17,8.33703e-18,-9.06646e-18,-3.8515e-18,-5.81143e-17,5.69069e-17,-1.82069e-17,-7.4988e-33,-6.18874e-18,-3.01082e-17,-2.21757e-17,4.39318e-17,5.04465e-17,-3.1586e-17,1.78732e-17,-4.06539e-17,-5.79085e-17,4.69064e-17,2.55971e-17,-1.39842e-16,5.47099e-17,-4.89026e-17,8.68092e-17,-5.11143e-17,-1.4769e-16,-3.12377e-16,3.36172e-16,1.8074e-16,9.37955e-17,-7.03973e-17,-2.28675e-16,-4.09591e-16,2.7095e-16,9.03613e-17,-6.1885e-17,-2.23909e-16,-3.75354e-16,3.30025e-16,1.67007e-16,6.0668e-21,-1.67817e-16,-3.25223e-16,3.77949e-16,2.23421e-16,5.74462e-17,-1.27365e-16,-2.58099e-16,4.15952e-16,2.27728e-16,1.16854e-16,-7.15707e-19,-2.0833e-16,-2.63123e-16,2.63309e-16,8.68445e-17,-2.21658e-17,-6.2127e-17,-1.62114e-16,2.75169e-16,1.21911e-16,5.22588e-17,-9.57627e-17,-1.18077e-16,-1.51817e-16,1.05898e-16,2.72126e-17,3.8823e-17,-8.32791e-17,-3.49387e-17,6.61191e-17,2.42159e-17,2.99952e-32,1.79688e-19,1.81045e-17,-1.06816e-16,-6.78235e-17,3.08193e-17,-3.96358e-18,5.1892e-17,-1.17234e-16,-1.05474e-16,-1.04636e-16,2.24657e-17,7.18047e-17,2.05405e-16,-2.50833e-16,-5.80291e-17,-3.54813e-17,1.67972e-16,2.1771e-16,-3.23515e-16,-1.71544e-16,-8.78146e-17,7.34154e-17,2.2899e-16,4.07471e-16,-2.75233e-16,-9.65304e-17,5.41127e-17,2.1482e-16,3.92993e-16,-3.27e-16,-1.64424e-16,-1.01113e-20,1.63461e-16,3.28248e-16,5.00112e-16,6.55668e-16,7.95204e-16,-7.66982e-16,-5.78852e-16,-4.18072e-16,-2.27413e-16,-1.13836e-16,6.69655e-18,2.17527e-16,2.75781e-16,4.1918e-16,4.88549e-16,6.01704e-16,-4.64204e-16,-3.59666e-16,-2.37099e-16,-7.89264e-17,-4.19598e-18,3.80332e-17,6.57168e-17,2.15974e-16,1.85912e-16,2.70305e-16,-1.79602e-16,-1.62657e-16,-9.40529e-17,-7.81227e-17,-3.0225e-17,-6.74892e-32,5.82937e-18,-6.10089e-18,-9.72594e-17,-1.30329e-16,-2.23107e-16,2.61558e-16,2.1141e-16,1.64099e-16,1.57835e-16,5.1343e-17,-7.05286e-17,-1.14789e-16,-2.43475e-16,-3.3761e-16,-5.25862e-16,6.77188e-16,5.88902e-16,4.32068e-16,3.10857e-16,1.62347e-16,8.18338e-17,-7.64335e-17,-2.29306e-16,-4.05351e-16,-5.53151e-16,-7.29968e-16,8.14082e-16,6.54693e-16,5.05302e-16,3.37853e-16,1.6878e-16,1.41559e-20,-1.66044e-16,-3.31273e-16,-4.89995e-16,-6.46579e-16,-8.15188e-16,7.1764e-16,5.83136e-16,4.20193e-16,2.27097e-16,1.10818e-16,-1.26774e-17,-2.26724e-16,-2.88439e-16,-4.35536e-16,-5.08831e-16,-6.26131e-16,4.35424e-16,3.26334e-16,1.99028e-16,1.46964e-16,6.71554e-17,-9.13261e-17,-1.24379e-16,-1.69109e-16,-2.55677e-16,-2.34755e-16,2.09358e-16,7.55265e-17,1.12022e-16,-2.0896e-17,3.6234e-17,1.19981e-31,-1.18384e-17,-5.90275e-18,7.92901e-17,1.06438e-16,1.93351e-16,-1.86085e-16,-2.52667e-16,-2.10964e-16,-9.91727e-17,1.94982e-18,7.56919e-18,1.57774e-16,2.81545e-16,3.70942e-16,4.4362e-16,-6.52762e-16,-5.6862e-16,-4.15712e-16,-2.98199e-16,-1.5315e-16,-7.5853e-17,7.94516e-17,2.29621e-16,4.0323e-16,5.48867e-16,7.23799e-16,-7.94099e-16,-6.63782e-16,-4.87664e-16,-3.34827e-16,-1.66197e-16,-1.82004e-20,1.68627e-16,3.34298e-16,5.07634e-16,6.65245e-16,8.07415e-16,-7.23809e-16,-5.87419e-16,-4.22313e-16,-2.26782e-16,-1.078e-16,1.86582e-17,2.3592e-16,3.01097e-16,4.51892e-16,5.29113e-16,6.50557e-16,-5.17667e-16,-4.04025e-16,-2.71981e-16,-1.0398e-16,-1.90925e-17,3.35966e-17,7.20183e-17,1.22244e-16,2.1442e-16,1.99205e-16,-2.39115e-16,-9.94182e-17,-1.89692e-17,8.89231e-18,-4.22431e-17,-1.8747e-31,1.78475e-17,1.79064e-17,-6.13208e-17,-8.25461e-17,-1.63594e-16,-3.33477e-16,1.82902e-16,-5.19327e-16,1.51533e-16,-8.32399e-16,-5.56321e-17,9.09465e-16,-2.08593e-16,8.16972e-16,-4.724e-16,6.28335e-16,-7.8393e-16,3.99356e-16,-1.15775e-15,1.43954e-16,-1.48444e-15,-8.24697e-17,1.37989e-15,-4.0111e-16,1.12075e-15,-7.1763e-16,8.01871e-16,-1.07573e-15,4.97781e-16,-1.43068e-15,1.63614e-16,-1.77633e-15,-1.64271e-16,1.43903e-15,-4.97517e-16,1.09245e-15,-7.99643e-16,7.29978e-16,-1.07363e-15,4.24433e-16,-1.32785e-15,1.04782e-16,1.52967e-15,-1.34094e-16,1.12954e-15,-4.68248e-16,7.82873e-16,-6.74984e-16,4.88887e-16,-7.3953e-16,2.3391e-16,-8.27183e-16,8.2052e-17,8.01289e-16,-1.3068e-16,4.79733e-16,-1.73163e-16,2.80433e-16,-1.75218e-16,1.2331e-16,-1.85106e-16,3.11132e-18,-6.27702e-17,2.69957e-31,-2.38565e-17,8.11123e-17,-1.78693e-16,1.69677e-16,-1.99229e-16,2.97927e-16,-2.24159e-16,4.72462e-16,-9.28712e-17,7.74669e-16,-7.32736e-18,-8.6648e-16,2.46663e-16,-7.8364e-16,5.0118e-16,-6.03908e-16,8.04212e-16,-3.83e-16,1.05938e-15,-1.34757e-16,1.49042e-15,8.54878e-17,-1.37957e-15,3.9899e-16,-1.12503e-15,7.66972e-16,-8.09644e-16,1.0944e-15,-5.07897e-16,1.4337e-15,-1.6797e-16,-1.77638e-15,1.66854e-16,-1.43601e-15,4.874e-16,-1.10153e-15,8.19626e-16,-7.36148e-16,1.12486e-15,-4.26554e-16,1.32816e-15,-1.01764e-16,-1.52369e-15,1.43291e-16,-1.11688e-15,4.84604e-16,-7.62591e-16,5.88389e-16,-4.60107e-16,7.72861e-16,-1.9584e-16,8.70168e-16,-3.39891e-17,-7.47996e-16,7.83199e-17,-5.26598e-16,2.42928e-16,-3.15983e-16,2.56484e-16,-1.47202e-16,1.67137e-16,-1.5115e-17,5.67611e-17,-3.67441e-31,2.98656e-17,-6.91086e-17,1.96662e-16,-1.45785e-16,2.28985e-16,-2.62378e-16,1.54394e-16,-5.36619e-16,1.45232e-16,-8.27962e-16,-4.07355e-17,8.23495e-16,-2.84733e-16,7.50309e-16,-5.2996e-16,5.79482e-16,-8.24494e-16,4.77667e-16,-1.07204e-15,2.36583e-16,-1.4964e-15,-8.85059e-17,1.37926e-15,-3.96869e-16,1.07381e-15,-7.60803e-16,8.17416e-16,-1.08531e-15,4.90259e-16,-1.43673e-15,1.65387e-16,-1.77633e-15,-1.69437e-16,1.43298e-15,-5.05039e-16,1.08287e-15,-8.11854e-16,7.42317e-16,-1.12058e-15,4.28674e-16,-1.32848e-15,9.87456e-17,1.51771e-15,-1.52488e-16,1.10422e-15,-3.89937e-16,7.42309e-16,-6.12815e-16,4.31327e-16,-8.06193e-16,2.68792e-16,-9.13153e-16,-1.40738e-17,8.05726e-16,-1.36982e-16,4.62441e-16,-2.01671e-16,2.4051e-16,-2.26727e-16,1.71093e-16,-1.49168e-16,2.71186e-17,-5.0752e-17,4.79923e-31,-3.58747e-17,5.7105e-17,-2.14632e-16,1.21893e-16,-2.58742e-16,3.3785e-16,-1.95651e-16,4.89754e-16,-8.65697e-17,7.70233e-16,-2.22239e-17,-8.91533e-16,2.11781e-16,-8.28e-16,4.47718e-16,-6.66077e-16,7.33753e-16,-4.61311e-16,1.0847e-15,-2.27386e-16,1.50238e-15,9.1524e-17,-1.37894e-15,3.94749e-16,-1.07809e-15,7.54634e-16,-7.97433e-16,1.07622e-15,-5.00376e-16,1.43975e-15,-1.69743e-16,-1.77639e-15,1.65081e-16,-1.42996e-15,4.94922e-16,-1.09196e-15,8.04082e-16,-7.48486e-16,1.11629e-15,-4.30794e-16,1.32879e-15,-9.57275e-17,-1.51173e-15,1.61684e-16,-1.09156e-15,4.06293e-16,-8.33049e-16,6.37242e-16,-5.13569e-16,7.28502e-16,-2.30722e-16,8.45115e-16,-4.88856e-17,-7.52433e-16,8.46214e-17,-5.09306e-16,1.60414e-16,-2.7606e-16,1.96971e-16,-8.39626e-17,2.42221e-16,-3.91222e-17,4.4743e-17,-6.07403e-31,4.18837e-17,-4.51014e-17,2.32601e-16,-9.80017e-17,1.77476e-16,-3.02301e-16,2.36908e-16,-5.53911e-16,1.3893e-16,-8.23526e-16,-2.5839e-17,8.48549e-16,-2.49852e-16,7.94668e-16,-4.76498e-16,6.41651e-16,-7.54035e-16,4.44955e-16,-1.09736e-15,2.18189e-16,-1.50836e-15,-9.45421e-17,1.37862e-15,-3.92629e-16,1.08237e-15,-7.48464e-16,8.05205e-16,-1.09488e-15,4.82737e-16,-1.44278e-15,1.6716e-16,-1.77632e-15,-1.67664e-16,1.44081e-15,-4.84805e-16,1.10105e-15,-7.9631e-16,7.54655e-16,-1.11201e-15,4.32915e-16,-1.32911e-15,9.27094e-17,1.50575e-15,-1.70881e-16,1.0789e-15,-4.22649e-16,8.12767e-16,-6.61669e-16,4.84789e-16,-7.61833e-16,1.92652e-16,-8.88099e-16,8.22762e-19,8.10162e-16,-1.43283e-16,5.56171e-16,-2.30179e-16,3.1161e-16,-2.78237e-16,1.07854e-16,-2.24251e-16,5.11259e-17,-3.87339e-17,7.4988e-31,-4.78928e-17,3.30977e-17,-1.39548e-16,1.85132e-16,-2.07233e-16,2.66751e-16,-1.67143e-16,5.07045e-16,-8.02682e-17,7.65796e-16,7.39018e-17,1.08182e-15,-1.93252e-15,-7.61337e-16,5.05278e-16,1.93629e-15,-1.89022e-15,-4.28599e-16,1.11001e-15,2.78861e-15,-1.59428e-15,9.75601e-17,1.84134e-15,-2.88465e-15,-1.08666e-15,7.42295e-16,2.62871e-15,-2.38365e-15,-4.92854e-16,1.43192e-15,3.3812e-15,-1.7764e-15,1.63308e-16,2.10105e-15,-3.02251e-15,-1.08238e-15,8.16293e-16,2.68087e-15,-2.22294e-15,-4.35035e-16,1.32942e-15,3.01893e-15,-1.49977e-15,1.80078e-16,1.70931e-15,-2.33655e-15,-7.92485e-16,5.75073e-16,1.87546e-15,-1.42528e-15,-2.65604e-16,8.20062e-16,1.8236e-15,-7.56869e-16,9.09229e-17,8.40254e-16,-1.03232e-15,-2.36137e-16,2.4848e-16,5.34388e-16,-3.4883e-16,-6.31295e-17,3.27249e-17,-9.07355e-31,5.39018e-17,-2.10941e-17,-2.86572e-16,6.15915e-16,2.36989e-16,-3.42224e-16,-1.01285e-15,8.72087e-16,1.32629e-16,-8.19089e-16,1.76541e-15,8.73602e-16,-2.1497e-16,-1.49244e-15,1.90843e-15,5.92797e-16,-7.94599e-16,-2.36331e-15,1.76391e-15,1.99796e-16,-1.52033e-15,3.00805e-15,1.37799e-15,-3.88388e-16,-2.23973e-15,2.65005e-15,8.20749e-16,-1.07671e-15,-3.02199e-15,2.10389e-15,1.68933e-16,-1.77631e-15,3.37988e-15,1.43476e-15,-4.92327e-16,-2.40573e-15,2.63317e-15,7.66993e-16,-1.10344e-15,-2.89351e-15,1.8344e-15,8.66732e-17,-1.61484e-15,2.80833e-15,1.16461e-15,-4.5536e-16,-1.89233e-15,1.84299e-15,5.38251e-16,-8.28497e-16,-1.88189e-15,1.13536e-15,1.57193e-17,-8.50736e-16,1.40473e-15,5.38879e-16,-2.58687e-16,-7.27514e-16,6.69455e-16,1.55638e-16,-1.88313e-16,-3.68956e-16,8.43065e-17,1.07983e-30,5.11114e-17,-3.23976e-16,-1.75487e-16,1.37349e-16,6.21433e-16,-6.92527e-16,-2.49657e-16,5.24337e-16,1.36932e-15,-9.03975e-16,5.90053e-17,1.05676e-15,-1.85638e-15,-8.05696e-16,4.51816e-16,1.87412e-15,-1.84965e-15,-3.95887e-16,1.13533e-15,2.807e-15,-1.58232e-15,1.03596e-16,1.84197e-15,-2.88889e-15,-1.09523e-15,7.29957e-16,2.64093e-15,-2.40183e-15,-4.85332e-16,1.43797e-15,3.37942e-15,-1.77641e-15,1.68474e-16,2.09322e-15,-3.01499e-15,-1.10056e-15,8.00748e-16,2.66853e-15,-2.23151e-15,-3.83764e-16,1.33005e-15,3.02497e-15,-1.48781e-15,1.98471e-16,1.73463e-15,-2.30384e-15,-7.51921e-16,6.23926e-16,1.93302e-15,-1.46964e-15,-1.89464e-16,9.06031e-16,1.8087e-15,-7.61306e-16,9.72244e-17,8.57545e-16,-1.00382e-15,-3.07236e-16,1.88967e-16,5.97627e-16,-3.84768e-16,-8.71368e-17,1.31729e-16,-1.2673e-30,6.59199e-17,2.91319e-18,-3.61656e-16,5.52677e-16,1.8548e-16,-2.71124e-16,-1.04135e-15,8.54795e-16,1.26327e-16,-8.14652e-16,1.78031e-15,8.98655e-16,-1.80088e-16,-1.44808e-15,1.85087e-15,6.54966e-16,-8.35163e-16,-2.285e-15,1.73859e-15,1.81403e-16,-1.42126e-15,3.05752e-15,1.37736e-15,-3.84147e-16,-2.23116e-15,2.66239e-15,8.08538e-16,-1.08628e-15,-3.02951e-15,2.09784e-15,1.63768e-16,-1.7763e-15,3.38166e-15,1.44259e-15,-4.99849e-16,-2.38755e-15,2.64872e-15,7.2382e-16,-1.09487e-15,-2.88927e-15,1.83377e-15,8.0637e-17,-1.51578e-15,2.78993e-15,1.1393e-15,-3.7705e-16,-1.82187e-15,1.90516e-15,4.80691e-16,-7.84137e-16,-1.84701e-15,1.16041e-15,3.06159e-17,-8.46299e-16,1.39843e-15,5.21587e-16,-1.76173e-16,-7.67437e-16,6.17945e-16,9.23987e-17,-1.52374e-16,-3.44949e-16,9.63246e-17,1.46976e-30,1.50116e-16,-3.47984e-16,-2.11425e-16,8.95656e-17,6.72942e-16,-7.63626e-16,-2.21149e-16,5.41629e-16,1.37562e-15,-9.08412e-16,4.41087e-17,1.14273e-15,-1.89127e-15,-7.39033e-16,5.09376e-16,1.92297e-15,-1.80909e-15,-4.74198e-16,1.16065e-15,2.71437e-15,-1.57036e-15,1.09633e-16,1.8426e-15,-2.89313e-15,-1.10379e-15,7.17619e-16,2.62538e-15,-2.39225e-15,-5.05566e-16,1.43015e-15,3.38459e-15,-1.77642e-15,1.66701e-16,2.09927e-15,-3.00747e-15,-1.09098e-15,8.12959e-16,2.65619e-15,-2.24008e-15,-3.88005e-16,1.33069e-15,3.03101e-15,-1.47585e-15,2.16864e-16,1.75994e-15,-2.27113e-15,-8.2238e-16,6.7278e-16,1.87956e-15,-1.514e-15,-2.24346e-16,8.80978e-16,1.7938e-15,-7.65742e-16,1.03526e-16,8.74837e-16,-9.75307e-16,-2.67313e-16,2.40476e-16,5.49843e-16,-3.09685e-16,-1.21761e-19,1.19711e-16,-1.68723e-30,7.79381e-17,-8.41018e-17,-3.25717e-16,6.0046e-16,2.44993e-16,-3.11047e-16,-9.58839e-16,8.37504e-16,1.20026e-16,-8.10216e-16,1.79521e-15,8.12686e-16,-2.56228e-16,-1.51474e-15,1.90433e-15,6.06113e-16,-7.64705e-16,-2.31772e-15,1.71328e-15,1.6301e-16,-1.43323e-15,3.05148e-15,1.37673e-15,-4.35418e-16,-2.22259e-15,2.67473e-15,7.96327e-16,-1.09586e-15,-3.00928e-15,2.10566e-15,1.65541e-16,-1.77629e-15,3.38343e-15,1.43654e-15,-5.07371e-16,-2.39713e-15,2.6365e-15,7.36158e-16,-1.08631e-15,-2.88503e-15,1.83313e-15,7.46008e-17,-1.52774e-15,2.77154e-15,1.11398e-15,-4.09761e-16,-1.86244e-15,1.85631e-15,5.34153e-16,-7.39778e-16,-1.92315e-15,1.07444e-15,4.55124e-17,-8.41863e-16,1.39212e-15,5.04295e-16,-2.04681e-16,-6.96338e-16,6.77458e-16,1.40182e-16,-2.27458e-16,-3.20942e-16,1.08343e-16,];
p.plot(X,Y) 
p.show()

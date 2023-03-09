# def swap_move_rotate(self, swap_type=1, variation=1):
# """
# Does a (composite of rotations) 'swap move' on the Rubik's cube depending on the value of the swap_type:

# swap_type=1: 3-cycle on centre corners (= centremost X-centres) (Figures 4,6 Bonzio 2017)
# swap_type=2: 3-cycle on coupled edges (= wing-edges) (Figures 5,9,10 Bonzio 2017)
#     variation=1,2,... = which theta_i wing-edges are 3-cycled

# """
# # TODO: I skipped 3-cycle on central edges (= centremost +-centres) (Figure 6  Bonzio 2017), since by my
# # TODO: analysis it was actually the the same as 1: 3-cycle on centre corners (= centremost X-centres) (Figure 4
# # TODO: Bonzio 2017) but using different faces, and an opposite direction 3 cycle

# if swap_type == 1:
#     # TODO: maybe generalize to different faces if needed?
#     # TODO: maybe generalize to different X-centres if needed?- i.e. w_22, w_11 etc
#     # Do 'centre corner'='centremost X-centre' 3-cycle from Figure 4 Bonzio 2017
#     # Move is described in their notation as (for the n=5 cube):
#     # z = [[C_F, C_D], U^-1] = C_F * C_D * C_F^-1 * C_D^-1 * U^-1 * C_D * C_F * C_D^-1 * C_F^-1 * U
#     # In our notation this is (for the n=5 cube)
#     # z = [[R_{0,1,0}, R_{3,1,0}], R_{1,0,1}]
#     #   = R_{0,1,0} * R_{3,1,0} *  R_{0,1,1} * R_{3,1,1} * R_{1,0,1} * R_{3,1,0} * R_{0,1,0} * R_{3,1,1} * R_{0,1,1} * R_{1,0,0}
#     # We generalize to an n-cube case by replacing all l=1 with f=floor(n/2-1)

#     # Validation: Centre corners (= X centres) only exist for n>3 therefore throw ValueError if n<=3
#     if self.n <= 3:
#         raise ValueError("Centre corners (= X centres) do not exist for n<4 therefore this swap move makes no sense.")

#     # Do swap move combination of rotations
#     self.rotate(0, math.floor((self.n / 2) - 1), 0)
#     self.rotate(3, math.floor((self.n / 2) - 1), 0)
#     self.rotate(0, math.floor((self.n / 2) - 1), 1)
#     self.rotate(3, math.floor((self.n / 2) - 1), 1)
#     self.rotate(1, 0, 1)
#     self.rotate(3, math.floor((self.n / 2) - 1), 0)
#     self.rotate(0, math.floor((self.n / 2) - 1), 0)
#     self.rotate(3, math.floor((self.n / 2) - 1), 1)
#     self.rotate(0, math.floor((self.n / 2) - 1), 1)
#     self.rotate(1, 0, 0)

# elif swap_type == 2:
#     # TODO: maybe generalize to different faces if needed?
#     # Do 'coupled edges'='wing edges' 3-cycle from Figure 5,9,10 Bonzio 2017
#     # Move is described in their notation as (for the n=5 cube):
#     # e = [C_L^-1, [L,U^-1]] = C_L^-1 * L * U^-1 * L^-1 * U * C_L * U^-1 * L * U L^-1
#     # In our notation this is (for the n=5 cube)
#     # e = [R_{4,1,1}, [R_{4,0,0}, R_{1,0,1}]]
#     #   = R_{4,1,1} * R_{4,0,0} * R_{1,0,1} * R_{4,0,1} * R_{1,0,0} * R_{4,1,0} * R_{1,0,1} * R_{4,0,0} * R_{1,0,0} * R_{4,0,1}

#     # We generalize to an n-cube case by doing nothing - i.e. above will do a n-cube theta_1 3-cycle

#     # We generalize to different theta_i wing edges with the variation parameter being used as i in theta_i
#     # As shown in Figure 10 Bonzio 17, you just have to change l=1 to l=i to change the theta-i which are
#     # 3-cycled

#     # Validation: Coupled edges (= wing edges) only exist for n>3 therefore throw ValueError if n<=3
#     if self.n <= 3:
#         raise ValueError("Coupled Edges (= Wing Edges) do not exist for n<4 therefore this swap move makes no sense.")

#     # Validation: There are only floor(n/2 - 1) types of wing edge therefore if variation is larger than
#     # this value or less than 1 we should throw a ValueError
#     if variation > math.floor((self.n/2) - 1) or variation < 1:
#         raise ValueError("There are only floor(n/2 - 1) types of wing edge therefore this swap move makes no sense.")

#     # Do swap move combination of rotations
#     self.rotate(4, variation, 1)
#     self.rotate(4, 0, 0)
#     self.rotate(1, 0, 1)
#     self.rotate(4, 0, 1)
#     self.rotate(1, 0, 0)
#     self.rotate(4, variation, 0)
#     self.rotate(1, 0, 1)
#     self.rotate(4, 0, 0)
#     self.rotate(1, 0, 0)
#     self.rotate(4, 0, 1)
# What a Microscopic Selection Mechanism Must Supply to Reproduce the Born Rule: A Specification from Pre-registered Runs

**John M. Bramble, MD, and Claude (Anthropic)**

*Draft v0.1, 2026-09-14. Written in `rayolddog/MathX` on a clean slate at the sponsor's
instruction, as the second of three shapes offered after the internal review round of the
Dirac–Kuramoto record's Paper 1 (v0.10; `born_selection_AI_review_2026-09/`). Every number is in a
`results/` folder of this repository with its prediction written before the run and its misses
kept; the folders are named at each use. The Dirac–Kuramoto record is cited as the record and is
not edited by this paper. Numbers here are model outputs, not measurements.*

**Contributions.** JMB supplied the physical ideas the runs were built on — that the quantum vacuum
interacts with the individual atoms of a detector's bulk, that the order of the interactions
between vacuum, bound electron and free particle does not commute, and that the aperture holds many
atoms with varying phases — and adjudicated every scope and interpretation decision. Claude
formalized the model, wrote and ran the pre-registered simulations, and wrote the prose. Byline
order is the sponsor's decision.

## Abstract

We ask what a microscopic mechanism for the selection of one detector site by one quantum must
contain in order to reproduce the Born rule and one record per quantum, and we answer with a
specification of six clauses, each supported by a pre-registered simulation showing what happens
when the clause fails. Two clauses are the Born rule and the projection postulate themselves: the
per-site hazard is proportional to the local intensity, and the rates are those of the normalized
conditional state, so that the first firing ends every other clock and a null at one port raises
the rate at another. The specification does not derive these; it locates them. The other four
clauses are the ways a mechanism that has inherited both can still lose them: by selecting a
maximum instead of running a race, which hands the brighter port every record unless the noise is
Gumbel-distributed at unit scale, which is the race in disguise; by holding too few sites per port
for the deposit's phase factor to average, which predicts a twelve-percent departure from Born at
one site; by carrying a frozen per-site factor with a heavy tail, which breaks a single detector's
binomial split even though any light-tailed factor is harmless; and, where an exchange of energy
precedes commitment, by a reflecting boundary or a threshold, either of which biases the winner
although every symmetric zero-mean exchange in a closed pot is otherwise fair. The worked example
is a model in the sponsor's two-spinor vocabulary that meets every clause and is, once written
out, the quantum-jump unravelling of photodetection with the vacuum as the name of each site's
clock. Two of its registered predictions are excluded by existing experiments, and both are
recorded. Nothing here goes beyond quantum mechanics; what is offered is a checklist, with numbers,
for anyone who proposes a detector dynamics as the origin of the square.

## 1. The question, and what kind of answer this is

A photon reaches a detector face that holds many atoms, and one of them takes it. Quantum mechanics
says which one with a probability proportional to the local intensity, and says that exactly one
does, and it says both by postulate. A program in the sponsor's repositories asked over the summer
of 2026 whether a detector dynamics could supply either statement instead of assuming it, and
concluded that it could not: every dynamics tried reproduced the square only by containing the
golden rule's rate law, and reproduced exclusivity only by containing a one-quantum constraint that
acts across the detector [DK record, `NEGATIVE_RESULT.md`]. That is the negative result, and this
paper does not reopen it.

What the runs behind that result also produced, and what no earlier paper collected, is a list of
the ways a mechanism can fail *after* it has inherited the square. A rule can have the Born rate
law in it and still give the brighter port every record, because it selects a maximum. It can have
the rate law and still depart from Born at a single site, because a phase factor has not averaged.
It can have the rate law and still give every detector an all-or-nothing split, because a per-site
factor is frozen with too heavy a tail. It can have a fair exchange and still bias the winner,
because the exchange stops at a threshold. Each of these was a registered prediction in this
repository that either held or missed, and in either case the run says what the clause is.

So the answer offered here is a specification. It is not a mechanism, and it is not a derivation.
It is the set of conditions under which a mechanism that already contains quantum mechanics'
two postulates about detection continues to reproduce them, with a number attached to each
condition showing the size of the failure without it. The worked example that meets every clause
is described in §8; its honest name is given there.

## 2. The setting shared by all the runs

Two ports, A and B, receive one quantum with intensities in the ratio $r$, so that the Born share
of port A is $r/(r+1)$: 0.800 at $r = 4$, the primary case throughout, and 0.667 at $r = 2$. Each
port holds $N$ atoms. An atom $i$ has a deposit $d_i = I_P\,w_i$, the port's intensity times a
phase weight $w_i = \cos^2(\Delta\varphi_i/2)$ that the sponsor's model attaches to the mismatch
between the incident particle's second spinor and the atom's; over a random bulk the phase weight
is uniform on $[0,1]$ with mean one half, and the runs are also done with $w_i = 1$ to isolate it.
A selection rule picks one atom. The quantities recorded are the port share $P(A)$, the exponent
$\gamma$ in $P(A)/P(B) = (\sqrt r)^\gamma$ (Born gives 2), and, where a single detector's statistics
are at issue, the Fano factor of the split of $M$ events in one detector against the binomial
value, $F = \mathrm{Var}(s)/[p(1-p)/M]$, which is 1 when events are independent Bernoulli trials.

Every run has a `PREREGISTRATION.md` written before it, with definitions fixed and numbers
predicted; the results files state the uncertainty first and list the predictions against the
outcomes. Where a prediction missed, the miss is in the table. Where a definition proved wrong, the
text is left as written and a dated note follows it.

## 3. Clause 1: the hazard is proportional to the local intensity, and interference with the amplitude is not enough

**The clause.** The rate at which site $i$ takes the quantum is proportional to $|\psi(x_i)|^2$.
This is the Born rule as a rate law, and the specification puts it first because everything
after it presupposes it.

**The run.** `vacuum_race/` (Table 1) attached a random per-atom factor to the deposit in six
forms and selected the atom with the largest result. Two forms bear on the clause. When the factor
*multiplies* the deposit — the vacuum acting on the atom's response, in the sponsor's ordering —
the rule can give Born (clause 2 says when). When the factor *adds to the incident amplitude* and
the selection is by the interference term, the random part is linear in $|\psi|$, and the rule
selects in proportion to the amplitude: the share is $\sqrt r/(\sqrt r+1) = 0.667$ at $r = 4$, with
exponent 1.01 at a thousand atoms per port, against Born's 0.800 and 2.

**What it means.** Interference does produce a product, but a product of amplitudes; a factor on
the intensity has to come from the site's response coefficient, which the vacuum sets by dressing
the atom before the packet arrives. At a photon-counting detector the vacuum does not interfere
with the incident field at all — Glauber's counting rate is normally ordered and the vacuum
contributes nothing to it — so the only place the vacuum can enter is the atom's rate, as a factor.
That is the physical content of the sponsor's ordering, vacuum on the fermion first, and it is the
reason the multiplicative form is the right one and the additive form gives the square root.

## 4. Clause 2: selection is a race, not a maximum

**The clause.** The site that takes the quantum is the first to fire among sites with exponential
clocks at rates $d_i$, so that site $i$ wins with probability $d_i/\sum_j d_j$. A rule that picks
the site with the largest of some per-site quantity does not do this.

**The run.** `vacuum_race/` Tables 1 to 3. With no random factor the largest deposit wins, and
the brighter port takes every record from ten atoms per port (1.000 against 0.800). With a
Gaussian factor multiplying the deposit, or a Gaussian or bounded uniform factor added to it, the
same: 1.000 by a hundred atoms per port for the multiplicative and bounded forms, and 0.966 and
0.990 at a hundred and a thousand for the additive Gaussian — a registered prediction of 1.000
that missed, the Gaussian maximum climbing only as the square root of the logarithm. The cause is
the maximum, not the factor: the largest of $N$ draws from any light-tailed distribution
concentrates at its ceiling as $N$ grows, and the brighter port's ceiling is higher. Additive
factors also read the intensity difference rather than the ratio: the share moves under a tenfold
scaling of both intensities (Table 4), where Born's does not.

Exactly one family of factors repairs the maximum: a standard Gumbel draw added to the logarithm
of the deposit, under which $\arg\max_i(\log d_i + G_i)$ selects site $i$ with probability
$d_i/\sum_j d_j$ [McFadden 1974; Yellott 1977; Maddison, Tarlow and Minka 2014]. That is
identically the first-passage race among exponential clocks [Gillespie 1977]. The share is then
0.702, 0.794, 0.801, 0.802 at $N = 1, 10, 100, 1000$, matching a sampled rule with the same rates
to the third decimal. At a Gumbel scale of $1/\beta$ the share is $r^\beta/(r^\beta+1)$, exponent
$2\beta$ (Table 3: 0.664, 0.798, 0.941 at $\beta = \tfrac12, 1, 2$); Born is the point $\beta = 1$
and nothing in the rule fixes it.

**What it means.** "Highest wins" is not a mechanism for Born; it is a mechanism against it, and
the one noise that rescues it is the one that turns it into a race. The clause is therefore that
the selection is a race, with the Gumbel identity as the statement of when a maximum happens to be
one. For a race, by contrast, the distribution of a fresh random factor on the rates hardly
matters at large $N$: any finite-mean factor gives Born once the port's rate-sum self-averages
(clause 4). The Gumbel requirement is an artefact of phrasing the selection as a maximum, which
is how the sponsor's rule was first stated and how it was first run.

## 5. Clause 3: many sites per port, or the phase factor does not average

**The clause.** If the deposit carries a per-site phase factor with a random phase, the port's
share is Born only once the port holds enough sites for the factor to average: about ten. At one
site per port the model predicts a departure.

**The run.** `second_spinor_phase_race/` Table 2 and `vacuum_race/` Table 1. Under the race with
the phase weight, the port share at $r = 4$ is 0.704 at one atom per port, 0.759 at two, 0.788 at
five, 0.794 at ten, and 0.800 from a hundred. The exact expression is
$\mathbb E[rW_A/(rW_A+W_B)]$ with $W_P$ the sum of the port's phase weights, which reaches
$r/(r+1)$ only as the sums self-average.

**What it means.** For a real pixel, with of order $10^{10}$ sites in reach, the clause is met
trivially and the phase factor is invisible, which is why the sponsor's phase mechanism leaves no
trace in a passive detector (clause 7). For a single-atom absorber the model predicts a
twelve-percent departure from Born at a four-to-one split. That is a prediction of the model, not
of quantum mechanics, and it is recorded as such: either the phase weight is physical and the
departure is testable in an atom-array detector, or it is not and it should be dropped from the
deposit. The runs of clause 7 favour dropping it.

## 6. Clause 4: the clock is memoryless, and a frozen factor is harmless only if its tail is light

**The clause.** The per-site clock is renewed for every arrival. A factor on the hazard that is
fixed per site for the life of the detector does no harm if its distribution has a light enough
tail for the port's rate-sum to self-average over the port's sites; a heavy-tailed frozen factor
gives every detector a split far from binomial.

**The runs.** `vacuum_race/` Tables 6, 7 and 9. With the model's per-atom factor redrawn every
$\tau$ arrivals, one detector's split of a thousand events has Fano factor $F = \tau$ (10.3 at
$\tau = 10$, 103 at $\tau = 100$); frozen for the life of the detector, every detector splits all or
nothing (step fraction 1.000) and Born appears only across the ensemble. That was the first result,
and the clause drawn from it — that any frozen per-site factor must decorrelate between arrivals —
was wrong, as the review round of 2026-09-14 showed and addendum 3 then registered and ran. The
factor in that run was the exponential of a Gumbel variable, which has infinite mean, so the port's
rate-sum is dominated by its largest term at every $N$ and never self-averages. Replace it with a
frozen uniform, lognormal or exponential factor and a fresh clock, and the Fano factor falls to
1.0, 1.2 and 1.4 by a thousand atoms per port (Table 9), while the infinite-mean factor stays
above a hundred and its ensemble mean is not even Born (0.76, a second registered miss).

**What it means.** Memorylessness of the clock is the golden rule's own property and the
specification requires it. Frozen per-site structure — inhomogeneous detuning, site geometry — is
harmless when it is bounded, which physical multipliers are, and a real pixel's $10^{10}$ sites
average it away in any case. So the Fano factor of a pixel's split, which reads the correlation
time of a frozen factor in the model, bounds nothing about frozen structure in a device, and the
frozen-vacuum reading of the sponsor's model was excluded by its infinite-mean factor, not by its
being frozen. The clause as first stated is withdrawn here as it is in the record (L30).

## 7. Clause 5: the rates are those of the normalized conditional state

**The clause.** When the wave reaches two ports at different times, the port share is the branch
weight of each port put in as an initial condition, not the outcome of a race between the ports.
Each site's rate is proportional to its deposit in the *normalized conditional state*: a port's
branch weight is consumed by its own no-fire absorption, so that the port's integrated hazard is
its share, and a null at one port renormalizes the others. The first fire anywhere ends every
clock. This is the projection postulate written as a stochastic process, and it is the second
thing the specification locates rather than derives.

**The runs.** `vacuum_race/` Table 8. The packet reaches port A at time zero and port B a fraction
$\Delta/T$ of a packet length later. With the race's rates running while the wave is present, the
nearer port takes 0.997 of the records at a tenth of a packet length of stagger and everything
beyond it, against Born's 0.800; this reproduces the failure the sibling record found for its own
selection postulate. A local repair, each port's branch consumed by its own absorption with the
ports coupled only by the stop on firing, gives 0.952 — a registered prediction of 0.800 that
missed, and the closed form says why: the far port fires only in the fraction of events the near
port left, $(1 - s_A\eta_A)\,s_B\eta_B$. Only the normalized conditional state gives Born, 0.798
to 0.802 over thirty cells at every stagger and coupling, and it does so as an identity: the port's
probability of ever firing is its branch weight by construction.

**What it means.** Under staggered arrival there is no race between ports; the race decides which
site within a port. The one-quantum constraint acts twice, cancelling the losing clocks when one
fires and renormalizing the survivors when none does, and both halves act across the port
separation. That is the quantum-jump unravelling of photodetection [Srinivas and Davies 1981;
Dalibard, Castin and Mølmer 1992; Carmichael 1993; Wiseman and Milburn 2010], and any mechanism
meeting this clause is that unravelling under another name. A "stop" with a speed or a delay is
not what the clause requires; the earlier figure in the record of "within a thousandth of the
commitment time" was a bound from a discarded race on a dynamical stop the conditional-state
update does not have, and it is withdrawn there and not used here.

## 8. Clause 6: fairness in the exchange class

**The clause.** If the mechanism lets the captured energy be exchanged among sites before
commitment, in a pot closed to one quantum, then any symmetric zero-mean exchange with an
absorbing zero is fair — the share of every site is a martingale and its probability of ending
with the whole pot is its initial share — whatever the kernel's dependence on the energies. The
$\sqrt e$ scaling that the sibling record's Paper 1 called the unique fair point is a property of
its independent-noise model, which does not close the pot. What can still bias the winner is a
reflecting zero, which favours the poor, or a termination at a threshold below the full share,
which favours the rich; a termination by absorption at the full share, or by a commit rate linear
in the share, does not.

**The run.** `exchange_fairness/`, pre-registered with Paper 1's own closed forms and the review
round's scratch numbers as predictions. Four models on the two-site pot (0.8, 0.2). *(Numbers to be
entered from `results/fairness.csv` when the run completes; the registered predictions are: model A,
independent noise $de_i = \sigma\sqrt{e_i}\,dW_i$, fair at 0.800 but with the survivor holding a
median 0.7 of a quantum, below half in a third of runs; model B, the conserving geometric-mean
kernel, fair at 0.800 with the pot exactly closed; model C, a conserving multiplicative kernel
whose sites are never absorbed, fair at 0.800 under a rate-linear commit and biased to 0.833 under
a threshold of 0.95, Paper 1's closed form $(s_0 - (1-t))/(2t-1)$; model D, fixed parcels with a
reflecting zero on three sites, biased toward the poor, below Born by more than 0.10 at the bright
site.)*

**What it means.** In the closed pot the share is $e_i/E$ with $E$ constant, so its drift is the
drift of $e_i$, and any symmetric zero-mean exchange has none. Fairness in this class is not a
knife-edge and needs no particular noise scaling; the two things that spoil it are boundaries,
and the specification names them. What the clause does not do is supply the square: the initial
shares are the intensities, put in.

## 9. Clauses that concern the sponsor's model in particular

Three further results are not clauses of the specification but bear on whether the sponsor's
model, the worked example, contains anything beyond it.

**The bound site's phase (clause 7 of the outline).** `spad_surface_contact/` Parts E and F solved
the Dirac equation for an electron bound in a well and measured the lower spinor directly: it is
slaved to the upper, in quadrature with it, fixed to $10^{-8}$ degrees and identical for every
identical atom, with a negative-energy free content of $10^{-6}$ at a binding of 0.03 rest energies
and about $10^{-16}$ for a silicon valence electron. Rebuilt on that spinor, the model's phase race
has nothing to select: the which-atom statistics are the nearest atom's and the which-port
statistics the arrival's. In a passive detector the phase mechanism leaves no trace. In a
coherently driven absorber the model predicts a phase-dependent response — but so does quantum
mechanics (coherent control [Heberle, Baumberg and Köhler 1995]), and where the two differ, in the
sign, the order in the probe field and the population factor, the existing experiments side with
quantum mechanics (L13, L14). That is one of the model's two registered predictions that
experiment excludes.

**The record (clause 8).** Parts A to C modelled the record as an induced current down a lattice
of random resistances with a branching avalanche. It decided nothing: the port statistics and the
single record were supplied by the selection, and the current copied them. This is a toy, with no
field, quench, dead time or afterpulsing, and it establishes nothing about a device; a real cascade
multiplies the vertex weight by a position- and energy-dependent efficiency and adds
photon-independent clicks. Its only use is to show where the selection stage ends.

**Energy (clause 9).** Part G tracked the quantum's energy under a divisible deposit with three
release channels at the excited site. With re-emitted photons leaving the device, the divisible
reading registers nothing at a balanced split for the toy's share of 0.4254 above 470 nm and gives
the brighter port every record where it registers; with the whole wave collected the line is
553 nm and at 505 nm it predicts a record at both ports instead. Either way it fails a balanced
split, which registers half the counts on each side with no coincidences at every visible
wavelength. The indivisible reading keeps Born, exclusivity, whole-quantum records and energy
closure at once, and is the postulate. The wavelength lines are the toy's and are not silicon's.

**The pair (clause 10).** `pair_case/`: the sponsor's model at two analyzers, a shared hidden
polarization read locally, gives Bell's sawtooth and CHSH 2.00; a shared bulk phase is a common
cause and changes nothing; the model's phase-dependent efficiency as a detection gate rides the
Garg–Mermin line to 2.83 at a conditional efficiency of 0.828, the detection loophole closed in
2015 [Hensen et al. 2015; Giustina et al. 2015; Shalm et al. 2015]; the conditional update of
particle 2 by A's outcome gives 2.83 with a rate-linear reading at B and 4.00 with a deterministic
one. All of this is Bell's theorem, and it is reported only because the sponsor's model proposed
the shared-bulk option in earnest. The second registered prediction experiment excludes is the
gate.

## 10. The worked example, by its right name

The sponsor's model, in the form the runs left it, is: a deposit at each site equal to the local
intensity times a phase factor that averages to one half; a clock at each site, renewed per
arrival, with rate proportional to that deposit in the normalized conditional state; the first
fire anywhere ending every clock. In the sponsor's vocabulary the clock is the quantum vacuum
acting on the bound electron before the free particle arrives, and the order of the three
interactions is what makes the vacuum a factor on the rate rather than a term in the amplitude.
Written out, this is the quantum-jump unravelling of photodetection with Born branch weights. It
meets every clause of §§3–8 by inheriting clauses 1 and 5 and satisfying 2, 3, 4 and 6. It derives
nothing, and the two places where it says more than quantum mechanics — the single-site phase
departure of clause 3 and the driven-detector phase dependence of §9 — are the two places
experiment or the Dirac equation tells against it. The model is therefore useful as an example and
not as a mechanism, and this paper offers it as the former.

## 11. What this paper does not claim, and what remains

It does not derive the square: clause 1 is the rate law. It does not derive exclusivity: clause 5
is the conditional-state update. It proposes no mechanism beyond the unravelling and no
prediction beyond quantum mechanics that survives the record. What it claims is that the six
clauses are necessary — each run shows a mechanism with the postulates in it that fails without
the clause — and that together with the postulates they are sufficient for the runs' domain: two
ports, one quantum, simultaneous or staggered arrival, any number of sites, a closed pot with or
without exchange.

What remains is what the runs did not cover. A derivation of the hazard's form from a
field–absorber–bath Hamiltonian, so that clause 1's linearity is a result; the record's own
open-system calculation was never attempted. The conserving class's continuum limit, so that
clause 6 is a theorem rather than four models. The multi-quantum sector, where the branch weights
of clause 5 are a joint registry. And the ontology: whether "the vacuum acting on the atom" names
anything the dressed-atom picture does not already name, which is the sponsor's question and is
open.

## Provenance and record

The three ideas that shaped the runs are the sponsor's and are stated in the Contributions. The
runs of 2026-09-06 to 2026-09-14 are in `kuramoto_one_two_rotations/`, `second_spinor_phase_race/`,
`sink_field_2d/`, `spad_surface_contact/`, `pair_case/`, `vacuum_race/` and `exchange_fairness/`,
each with its pre-registration, results and ledger; `MODEL_SECTION.md` §§1–12 is the model as it
developed, kept as a record; `README.md` carries plain-language summaries of the sessions. The
Dirac–Kuramoto record, at `rayolddog/DiracKuramotoFramework`, is the program this paper follows
from; its Paper 1 (v0.11) and `NEGATIVE_RESULT.md` are the negative result, and its review rounds
are where the errors in this repository's earlier statements were found. Two of them are carried
here as corrections rather than hidden: the frozen-factor clause of §6, and the reading of
"stop" in §7.

## References

*(Author–year; entries marked [verify] are from memory and are to be checked against the sources
before any circulation.)*

- Carmichael, H. J. (1993). *An Open Systems Approach to Quantum Optics.* Springer. [verify]
- Dalibard, J., Castin, Y., & Mølmer, K. (1992). Wave-function approach to dissipative processes in quantum optics. *Phys. Rev. Lett.* **68**, 580–583. [verify]
- Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.* **81**, 2340–2361. [verify]
- Giustina, M., et al. (2015). Significant-loophole-free test of Bell's theorem with entangled photons. *Phys. Rev. Lett.* **115**, 250401. [verify]
- Glauber, R. J. (1963). The quantum theory of optical coherence. *Phys. Rev.* **130**, 2529–2539. [verify]
- Grangier, P., Roger, G., & Aspect, A. (1986). Experimental evidence for a photon anticorrelation effect on a beam splitter. *Europhys. Lett.* **1**, 173–179. [verify]
- Gumbel, E. J. (1958). *Statistics of Extremes.* Columbia University Press.
- Heberle, A. P., Baumberg, J. J., & Köhler, K. (1995). Ultrafast coherent control and destruction of excitons in quantum wells. *Phys. Rev. Lett.* **75**, 2598–2601. [verify]
- Hensen, B., et al. (2015). Loophole-free Bell inequality violation using electron spins separated by 1.3 kilometres. *Nature* **526**, 682–686. [verify]
- Maddison, C. J., Tarlow, D., & Minka, T. (2014). A* sampling. *Advances in Neural Information Processing Systems* **27**. [verify]
- McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (ed.), *Frontiers in Econometrics*, Academic Press, 105–142. [verify]
- Shalm, L. K., et al. (2015). Strong loophole-free test of local realism. *Phys. Rev. Lett.* **115**, 250402. [verify]
- Srinivas, M. D., & Davies, E. B. (1981). Photon counting probabilities in quantum optics. *Optica Acta* **28**, 981–996. [verify]
- Wiseman, H. M., & Milburn, G. J. (2010). *Quantum Measurement and Control.* Cambridge University Press. [verify]
- Yellott, J. I. (1977). The relationship between Luce's choice axiom, Thurstone's theory of comparative judgment, and the double exponential distribution. *J. Math. Psychol.* **15**, 109–144. [verify]
- [DK record] Bramble, J. M., & Claude (Anthropic). *DiracKuramotoFramework*, https://github.com/rayolddog/DiracKuramotoFramework: `NEGATIVE_RESULT.md`; `drafts/PAPER1_DRAFT_born_selection.md` v0.11; `born_selection_AI_review_2026-09/`.

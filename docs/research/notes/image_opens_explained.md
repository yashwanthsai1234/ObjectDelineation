## Six words

- **Image**: one lung picture. It is one file on the disk.
- **Outline**: the shape that the model draws on an image.
- **Card**: one image and one question. The model reads one card. Then it draws one outline.
- **Stack**: all the cards for one step.
- **Open**: to read an image file and unpack it. The CPU opens an image before it puts the image on a card.
- **Round**: one part of the learning. One step has 16 rounds.

## Start small: 2 images

We want 2 outlines for each image.
The model draws one outline from one card.
So the CPU makes 4 cards: A, A, B, B.

```
 ┌───────── STACK ─────────┐
 │ card 1: image A  ◄── CPU opens A
 │ card 2: image A  ◄── CPU opens A again
 │ card 3: image B  ◄── CPU opens B
 │ card 4: image B  ◄── CPU opens B again
 └─────────────────────────┘
   4 opens. Only 2 are necessary.
```

Now there is a second problem.
The learning has 2 rounds.
The training program asks the CPU for a new stack in each round.
The CPU makes the same 4 cards again.
The model does not use the second stack.

```
 Stack 1:  A A B B   ──► the model uses it
 Stack 2:  A A B B   ──► nobody uses it
           ─────────
           8 opens. Only 2 are necessary.
```

## The parts

```
┌──────────┐  image files  ┌──────────┐   cards    ┌──────────┐
│   DISK   │ ────────────► │   CPU    │ ─────────► │   GPU    │
│ holds    │               │ opens    │            │ model    │
│ images   │               │ images,  │            │ reads    │
│          │               │ makes    │            │ cards,   │
│          │               │ cards    │            │ draws    │
│          │               │          │            │ outlines │
└──────────┘               └──────────┘            └──────────┘
```

## Flow: one real step

```
 START
   │
   ▼
 CPU makes stack 1:  256 cards, 256 opens  ──► the model uses it
 CPU makes stack 2:  256 cards, 256 opens  ──► nobody uses it
   ...
 CPU makes stack 16: 256 cards, 256 opens  ──► nobody uses it
   │       the GPU waits. This takes about 20 seconds (estimate).
   ▼
 GPU: the model reads 256 cards and draws 256 outlines
   │
   ▼
 GPU: the checker gives a score to each outline
   │
   ▼
 GPU: the model learns in 16 rounds
   │
   ▼
 END
```

## Interaction: who asks whom

```
 Training program        CPU                Disk          GPU (model)
       │                  │                  │                 │
       │── stack 1? ─────►│── open 256 ─────►│                 │
       │◄── 256 cards ────│                  │                 │
       │── stack 2? ─────►│── open 256 ─────►│                 │
       │◄── 256 cards ────│  (nobody uses)   │                 │
       │    ... stacks 3 to 16, the same ... │                 │
       │                                                       │
       │── 256 cards from stack 1 ────────────────────────────►│
       │◄── 256 outlines ──────────────────────────────────────│
```

## Example 1: the small case

```
 Images:                  2
 Outlines for each image: 2
 Cards in one stack:      2 × 2 = 4
 Stacks:                  2
 Opens:                   4 × 2 = 8
 Opens necessary:         2
```

## Example 2: the real step

```
 Images:                  32
 Outlines for each image:  8
 Cards in one stack:      32 × 8  = 256
 Stacks:                  16
 Opens:                   256 × 16 = 4,096
 Opens necessary:         32
```

## Example 3: one image

```
 Cards with this image in one stack:   8
 Stacks:                               16
 Opens of this image:                  8 × 16 = 128
 Opens necessary:                      1
```

We counted the opens on this Mac with the real program and 200 real images. The count was 4,096.

## Summary

The model draws one outline from one card, so the CPU makes 8 cards for each image.
The CPU opens the image again for each card, and it makes the same stack 16 times.
Result: 4,096 opens in each step. Only 32 are necessary.

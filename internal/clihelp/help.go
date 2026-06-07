// Package clihelp provides small helpers for grouped command-line help screens.
package clihelp

import (
	"fmt"
	"io"
	"strings"
)

// Flag describes one logical help item. Names may include aliases that are
// parsed as separate flags but should be documented together.
type Flag struct {
	Names   []string
	Arg     string
	Default string
	Text    string
}

// Group describes a titled block of related help items.
type Group struct {
	Title string
	Intro []string
	Items []Flag
}

// WriteFlagGroups emits consistently formatted grouped flag documentation.
func WriteFlagGroups(w io.Writer, groups []Group) {
	for _, group := range groups {
		_, _ = fmt.Fprintf(w, "%s:\n", group.Title)
		for _, line := range group.Intro {
			_, _ = fmt.Fprintf(w, "  %s\n", line)
		}
		if len(group.Intro) > 0 && len(group.Items) > 0 {
			_, _ = fmt.Fprintln(w)
		}
		for _, item := range group.Items {
			writeFlag(w, item)
		}
		_, _ = fmt.Fprintln(w)
	}
}

// WriteSizeModelReference emits the model-specific size-selection argument
// reference shared by radigest and radigest-design.
func WriteSizeModelReference(w io.Writer, dash string) {
	if dash == "" {
		dash = "-"
	}
	_, _ = fmt.Fprintln(w, "Size-selection models:")
	_, _ = fmt.Fprintln(w, "  hard")
	_, _ = fmt.Fprintf(w, "      Uses only the hard %smin/%smax window. Extra model parameters: none.\n", dash, dash)
	_, _ = fmt.Fprintf(w, "      Example: %smin 300 %smax 600 %ssize-model hard\n", dash, dash, dash)
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "  normal")
	_, _ = fmt.Fprintf(w, "      Applies a Gaussian-like recovery weight over %sscore-min/%sscore-max.\n", dash, dash)
	_, _ = fmt.Fprintf(w, "      Uses %ssize-mean and %ssize-sd.\n", dash, dash)
	_, _ = fmt.Fprintf(w, "      Example: %smin 300 %smax 600 %sscore-min 1 %sscore-max 2000 %ssize-model normal %ssize-mean 275 %ssize-sd 85\n", dash, dash, dash, dash, dash, dash, dash)
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "  triangular")
	_, _ = fmt.Fprintf(w, "      Applies a triangular recovery weight. Uses %ssize-mean, which must be inside (%smin,%smax).\n", dash, dash, dash)
	_, _ = fmt.Fprintf(w, "      Example: %smin 300 %smax 600 %ssize-model triangular %ssize-mean 450\n", dash, dash, dash, dash)
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "  soft-window")
	_, _ = fmt.Fprintf(w, "      Applies soft logistic edges around %smin/%smax. Uses %ssize-edge-sd.\n", dash, dash, dash)
	_, _ = fmt.Fprintf(w, "      Example: %smin 300 %smax 600 %sscore-min 1 %sscore-max 2000 %ssize-model soft-window %ssize-edge-sd 25\n", dash, dash, dash, dash, dash, dash)
	_, _ = fmt.Fprintln(w)
}

func writeFlag(w io.Writer, item Flag) {
	names := strings.Join(item.Names, ", ")
	if item.Arg != "" {
		names += " " + item.Arg
	}
	if item.Default != "" {
		_, _ = fmt.Fprintf(w, "  %-42s default: %s\n", names, item.Default)
	} else {
		_, _ = fmt.Fprintf(w, "  %s\n", names)
	}
	for _, line := range wrapHelp(item.Text, 78) {
		_, _ = fmt.Fprintf(w, "      %s\n", line)
	}
}

func wrapHelp(text string, width int) []string {
	text = strings.TrimSpace(text)
	if text == "" {
		return nil
	}
	words := strings.Fields(text)
	if len(words) == 0 {
		return nil
	}
	lines := make([]string, 0, 1)
	var line strings.Builder
	for _, word := range words {
		if line.Len() == 0 {
			line.WriteString(word)
			continue
		}
		if line.Len()+1+len(word) > width {
			lines = append(lines, line.String())
			line.Reset()
			line.WriteString(word)
			continue
		}
		line.WriteByte(' ')
		line.WriteString(word)
	}
	if line.Len() > 0 {
		lines = append(lines, line.String())
	}
	return lines
}

Attribute VB_Name = "Module1"
Sub extracttxt()
For i = 1 To 452
Cells(i, 4) = Mid(Cells(i, 1), 70, 28)
Next
End Sub

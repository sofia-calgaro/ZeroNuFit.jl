module try_me

greet() = print("Hello World!")

end # module YourPackageName

using .try_me
try_me.greet()
Here are some notes:

1. In many contact cases bodies aren't interesced, so stiffness matrix of the contact part is zero, in this regard JF modification for Newton-Krylov is not working properly, and therefore, removed in contacts.
2. Newton-Krylov methods converge fast, but require lots of sub-steps. Additionally, they aren't exact, which makes them especially unhelpful for contact problems.
4.   
